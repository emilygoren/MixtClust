
                                        # Row match: return row numbers of matrix that have rows equal to row.
row_match <- function(row, matrix) {
    ans <- apply(matrix, 1, function(r) all.equal(row, r) == TRUE)
    return(which(ans == TRUE))
}

## Perform one EM iteration (of all ECM steps).
EM_iter <- function(oldpars, x, A, Ru, miss.grp, ps, sigma.constr, df.constr, approx.df, marginalization, Q2old = NULL) {
    k <- length(oldpars$pi)
    ## 1st cycle -- update cluster proportions, means, dfs
    ## E-step followed by CM-step for pi
    ##    print(are.sigmas.valid(oldpars$Sigma))
    w <- up_W(x, oldpars$mu, oldpars$Sigma, oldpars$nu, miss.grp, Ru)
    if (k > 1) {
        z <- up_Z(x, oldpars$mu, oldpars$Sigma, oldpars$nu, oldpars$pi, miss.grp, Ru)
        pisnew <- as.numeric(up_pi(z))
    } else {
        z <- matrix(1, nrow(x), k)
        pisnew <- 1
    }
    if (marginalization) {
        ## calculate Q2 before
        if (is.null(Q2old))
            Q2old <- Q2(x, z, w, oldpars$Sigma, oldpars$mu, miss.grp, Ru)
        ## update locations 
        musnew <- up_mu(x, z, w, A)
        Q2newMu <- Q2(x, z, w, oldpars$Sigma, musnew, miss.grp, Ru)
        if (Q2old > Q2newMu) {
            musnew <- oldpars$mu
            ##    cat("Q2old = \n", Q2old, "Q2newMu" = Q2newMu, "musnew no good\n")
        }
    } else {
        K <- ncol(z)
        M <- length(unique(miss.grp))
        p <- ncol(x)
        ## EM for missing components
        SOiOEOO <- lapply(1:K, function(k) SOiOEOOk(oldpars$Sigma[,,k], Ru))
        xhat <- lapply(1:K, function(k) xhatk(x, oldpars$mu[k,], miss.grp, M, SOiOEOO[[k]]))
                                        # update locations
        musnew <- up_mu_Lin(p, z, w, xhat)
    }
    ## update nu's
    nusnew <- as.numeric(up_nu(z, w, oldpars$nu, ps, df.constr, approx.df))
    bidx <- !is.finite(nusnew); nusnew[bidx] <- oldpars$nu[bidx] # fix any NaN
    ## 2nd cycle -- update dispersions
    if (k > 1)
        z <- up_Z(x, musnew, oldpars$Sigma, nusnew, pisnew, miss.grp, Ru)
    w <- up_W(x, musnew, oldpars$Sigma, nusnew, miss.grp, Ru)
    if (marginalization) {
        ## calculate Q2 before
        Q2old <- Q2(x, z, w, oldpars$Sigma, musnew, miss.grp, Ru)
        Sigmasnew <- up_Sigma(x, z, w, musnew, A, sigma.constr)
        Q2newSigma <- Q2(x, z, w, Sigmasnew, musnew, miss.grp, Ru)
        if (Q2old > Q2newSigma) {
            Sigmasnew <- oldpars$Sigma
            ## cat("Q2old = \n", Q2newMu, "Q2newSigma" = Q2newSigma, "Sigmasnew no good\n")
        }
    } else {
        Sigmasnew <- sapply(1:K, 
                            function(k) up_Sigmak_Lin(M, z[,k], w[,k], musnew[k,], oldpars$Sigma[,,k], xhat[[k]], miss.grp, SOiOEOO[[k]]), 
                            simplify = 'array')
        if (sigma.constr) {
            wtdSigmas <- lapply(1:K, function(k) pisnew[k]*Sigmasnew[,,k])
            Ss <- Reduce('+', wtdSigmas)
            for (k in 1:K) Sigmasnew[,,k] <- Ss
        }
    }
    ## Output
    newpars <- list(pisnew, nusnew, musnew, Sigmasnew)
    names(newpars) <- c('pi', 'nu', 'mu', 'Sigma')
    return(newpars)
}

## Function to generate a set of initial parameter values.
get.init.val <- function(X, R, K, df.constr, sigma.constr, init = "smart-random", Z = NULL, eps = 1e-3) {
                                        # Draw transformation parameter, degrees of freedom, proportions uniformly.
                                        #if (df.constr) nus <- rep(runif(1, 5, 25), K) else nus <- runif(K, 5, 25)
                                        # Follow teigen: set dfstart to 50
                                        #nus <- runif(K, min = 10, max = 50)
    nus <- rep(50, K)
    if (df.constr) nus <- rep(nus[1], K)
    n <- nrow(X); p <- ncol(X)

    if (init == "smart-random")
        old.inits <- FALSE else
                              old.inits  <- TRUE
        
    if (!old.inits) {    
        ## smart random initialization chooses random points in sequence from the dataset with probability inversely proportional to the distance from (the closest of) the previous points. We use the kmmeans with 0 iterations to achieve this. Note that this function can handle missing data and so does not particularly care over missing cases, and can handle them fine.
        ## because we check so much, there is a chance we may never come out so we limit it to ten tries, otherwise we bail out and go over to EMG's uniform initialization
        y <- X
        y[R] <- NA
        Sigmas <- array(0, dim=c(p, p, K))
        if (K == 1) {
            Sigmas[,,1]  <- cov(y, use = "pairwise.complete.obs")
            pis <- 1
            mus <- matrix(colMeans(y,na.rm=T), nrow = 1)
        } else {
            sigmas.stab <- F
            grand.iter  <- 1
            while ((!sigmas.stab) & (grand.iter < 10)) {
                res <- kmmeans::kmmeans(data = as.data.frame(y), K = K, n.init = 1, kmmns.iter = 0)
                res$partition <- res$partition + 1
                nks <- table(res$partition)
                if (min(nks) > p) {
                    ## don't even bother otherwise (perhaps can and should be made even more stringent)
                    grand.iter  <- grand.iter  + 1
                    pis <- table(res$partition)/n
                    mus <- res$centers
                    for (k in 1:K) {
                        Sigmas[,,k] <- cov(y[res$partition==k,], use = "pairwise.complete.obs")
                        if (any(is.na(Sigmas[,,k])) | max(diag(Sigmas[,,k]) < eps*0.01))
                            break
                        eig.s <- eigen(Sigmas[,,k], symmetric=TRUE, only.values=TRUE)
                        if (eig.s$values[p] < (eps * eig.s$values[1])) {
                            Sigmas[,,k]  <- Sigmas[,,k] + diag(p)*eps*max(diag(Sigmas[,,k]))
                            eig.s <- eigen(Sigmas[,,k], symmetric=TRUE, only.values=TRUE)
                            if (eig.s$values[p] < (eps * eig.s$values[1]))
                                break ##not fixable easily, so bail out, try afresh
                        } else {
                            if (k == K)
                                sigmas.stab <- T
                        }
                    }
                }
            }
            if (!sigmas.stab)
                old.inits  <- TRUE
        }
        if (old.inits)
            init  <- "uniform"
    }
    if (old.inits) {
        minstable <- p/n
        CC <- rowSums(R) == 0
            ## Uniform initialization of z.
            if (K > 1 & is.null(Z)) {
                nonstable <- TRUE
                while(nonstable) {
                    Z <- switch(init,
                                uniform = {
                                    matrix(runif(n*K, 0, 1), n, K)
                                },
                                kmeans = {
                                    id <- kmeans(X[CC, ], nstart = 1e3, centers = K)
                                    sapply(1:K, function(k) as.numeric(id$cluster == k))
                                })
                    Z <- Z/rowSums(Z)
                    pis <- colMeans(Z)
                    if (all(pis > minstable)) nonstable <- FALSE
                }
            } else if (K > 1 & !is.null(Z)) {
                Z <- Z/rowSums(Z)
                pis <- colMeans(Z)
            } else {
                Z <- matrix(1, n, K)
                pis <- 1
            }
            ## Use updated IDs to set remaining parameters.
            mus <- matrix(0, K, p)
            Sigmas <- array(0, dim=c(p, p, K))
            for (k in 1:K) {
                wtk <- switch(init,
                              uniform = {
                                  Z[,k] * CC # set wt to zero if not complete case
                              },
                              kmeans = {
                                  out <- rep(0, n)
                                  out[CC] <- Z[,k]
                                  out
                              },
                              ids = {
                                  Z[,k] * CC
                              })
                wtcov <- cov.wt(X, wt = wtk, method = "ML")
                mus[k,] <- wtcov$center
                Sigmas[,,k] <- wtcov$cov
            }
    }
    if ((K > 1) & sigma.constr) {
        wtdSigmas <- lapply(1:K, function(k) pis[k]*Sigmas[,,k])
        S <- Reduce('+', wtdSigmas)
        for (k in 1:K) Sigmas[,,k] <- S
    }
    out <- list(pi = pis, nu = nus, mu = mus, Sigma = Sigmas)
    return(out)
}

## Run EM.
run.EM <- function(init, nclusters, X, miss.grp, A, Ru, ps, max.iter, tol, convergence, sigma.constr, df.constr, approx.df, marginalization, npar) {
    ##print("entering run.EM")
    old <- init
    del <- 1e6 # Initialize convergence holder.
    iter <- 0
    initL <- sapply(1:nclusters, function(k) {old$pi[k] * h(X, old$mu[k,], old$Sigma[,,k], old$nu[k], miss.grp, Ru)})
    LLs <- c(sum(log(rowSums(initL))), rep(NA, max.iter)) # Store loglikelihood at each iteration.
    while (del > tol) {
        iter <- iter + 1
        new <- EM_iter(old, X, A, Ru, miss.grp, ps, sigma.constr, df.constr, approx.df, marginalization)
        newL <- sapply(1:nclusters, function(k) {new$pi[k] * h(X, new$mu[k,], new$Sigma[,,k], new$nu[k], miss.grp, Ru)})
        newLLn <- sum(log(rowSums(newL)))
        ##    cat("ITER =", iter, "newLLn" = newLLn, "\n")
        LLs[iter+1] <- newLLn
        old <- new
                                        # Stop if loglik didn't change or went down
        if (newLLn <= LLs[iter]) break
                                        # Otherwise check specificed convergence critrion
        if (convergence == 'lop') {
            del <- ((LLs[iter] - newLLn) / LLs[iter]) # Relative change in loglikelihoods.
        } else if (convergence == 'aitkens') {
            if (iter > 1 ) {
                a <- (newLLn - LLs[iter]) / (LLs[iter] - LLs[iter - 1])
                LLinf <- LLs[iter] + ((newLLn - LLs[iter]) / (1 - a))
                del <- abs(LLinf - newLLn)
            } else {
                del <- 1e6
            }
        } else {
            stop('Specified convergence criterion not recognized.')
        }
        if (iter == max.iter) break
    }
                                        # Final posterior probabilities of cluster membership
    Zs <- up_Z(X, new$mu, new$Sigma, new$nu, new$pi, miss.grp, Ru)
    BIC <- 2*LLs[iter+1] - npar*log(nrow(X))
    res <- list(estimates = new, 
                iterations = iter, 
                Zs = Zs, 
                loglik = LLs[2:(iter+1)], 
                bic = BIC)
    ##    print("exiting run.EM")
    return(res)
}

## Run an initialization with short em -- don't keep track of LL, BIC, etc to save time
run.em <- function(nclusters, X, miss.grp, A, R, Ru, ps, niter, sigma.constr, df.constr, marginalization, init = "smart-random") {
    ## print("entering run.em")
    llhd = -Inf
    while (!is.finite(llhd))    {
        old <- get.init.val(X, R, nclusters, df.constr, sigma.constr, init)
        iter <- 0
        while (iter < niter) {
            new <- EM_iter(old, X, A, Ru, miss.grp, ps, sigma.constr, df.constr, approx.df = TRUE, marginalization)
            old <- new
            test <- are.sigmas.valid(new$Sigma)
            ## cat(iter,"in run.em", test, "\n")
            if (!test) {
                llhd <- -Inf
                iter <- niter
            } else {
                L <- sapply(1:nclusters, function(k) {new$pi[k] * h(X, new$mu[k,], new$Sigma[,,k], new$nu[k], miss.grp, Ru)})
                llhd <-  sum(log(rowSums(L)))
                iter <- iter + 1
            }
        }
    }
    ## Final loglik
    #    if (are.sigmas.valid(new$Sigma)) {
#    L <- sapply(1:nclusters, function(k) {new$pi[k] * h.try(X, new$mu[k,], new$Sigma[,,k], new$nu[k], miss.grp, Ru)})
#    llhd <-  sum(log(rowSums(L)))
    ##    } else
    ##       llhd <- -Inf
#}
    res <- list(estimates = new, loglik = llhd)
    ## print("exiting run.em")
    return(res)
}

##
## check for positive definiteness of the Sigmas
##

are.sigmas.valid <- function(Sigma, eps = 1e-4)
{
    all(sapply(1:dim(Sigma)[3], function(k) {
        if ((max(diag(Sigma[,,k])) > eps*0.01)  & (min(diag(Sigma[,,k]) > 0))) 
            eigen(cov2cor(Sigma[,,k]),symmetric = TRUE, only.values = TRUE)$values[dim(Sigma)[1]] > eps else FALSE
    }))
}

## Update Sigma -- R version used for debugging only
up_SigmaR <- function(x, z, w, mus, A, constr = FALSE) {
    k <- ncol(z)
    p <- ncol(x)
    n <- nrow(x)
    S <- sapply(1:k, function(kk) {
        muk <- mus[kk,]
        x.c <- sweep(x, 2, muk)
        Lk <- lapply(1:n, function(i) {
            Acent <- as.numeric(diag(A[i,]) %*% x.c[i,])
            return(z[i,kk] * w[i,kk] * Acent %o% Acent)
        })
        Lks <- Reduce("+", Lk)
        Rk <- lapply(1:n, function(i) {
            return((z[i,kk]) * A[i,] %o% A[i,])
        })
        Rks <- Reduce("+", Rk)
        return(Lks / Rks)
    }, simplify = 'array')
    if (constr) stop('Dispersion constraints not yet implemented')
    return(S)
}

                                        # Update mu -- R version used for debugging only
up_muR <- function(x, z, w, A) {
    k <- ncol(z)
    p <- ncol(x)
    n <- nrow(x)
    mu <- t(sapply(1:k, function(kk) {
        Wk <- lapply(1:n, function(i) z[i,kk] * w[i, kk] * diag(A[i,]))
        Lks <- solve(Reduce("+", Wk))
        Rk <- lapply(1:n, function(i) Wk[[i]] %*% x[i,])
        Rks <- Reduce("+", Rk)
        return(as.numeric(Lks %*% Rks))
    }))
    return(mu)
}

                                        # Loglik
loglikelihood <- function(x, ests) {
    X <- as.matrix(x)
    R <- is.na(X)
                                        # Remove observations with no observed coordinates and coordinates with no observations.
    bad.rows <- rowSums(R) == ncol(X)
    bad.cols <- colSums(R) == nrow(X)
    X[R] <- 0
    n <- nrow(X); p <- ncol(X)
                                        # Set missing values to 0 for cpp code.
                                        # Unique patterns of missingness.
    R.unique <- unique(R, MARGIN = 1)
                                        # Missingness pattern of each observation.
    miss.grp <- apply(R, 1, row_match, R.unique)
    Ru <- 1*!R.unique
    nclusters <- length(ests$pi)
    
    L <- sapply(1:nclusters, function(k) {ests$pi[k] * h(X, ests$mu[k,], ests$Sigma[,,k], ests$nu[k], miss.grp, Ru)})
    LL <- sum(log(rowSums(L)))
    return(LL)
}


#h.try <- function(x, mu, Sigma, nu, miss.grp, Ru) {
#    try.h <- try(h(x, mu, Sigma, nu, miss.grp, Ru))
#    print(try.h)
#    print(length(class(try.h))) # appears to be matrix(array)
#    print(length(class(try.h)))
#    if (length(class(try.h)==1)) {
#        if (class(try.h) == "try-error")
#            browser() 
#    }  else
#        try.h
#    h(x, mu, Sigma, nu, miss.grp, Ru)
#}

