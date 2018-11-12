
# use iris data w/ MCAR @ 10% missingness
d <- subset(iris, select = -Species)
set.seed(10)
R <- (matrix(rbinom(n = ncol(d)*nrow(d), size = 1, prob = 0.1), ncol = ncol(d)) == 1)
X <- as.matrix(d); X[R] <- NA
# clone X into Y and replace NA w/ zero
Y <- X; Y[R] <- 0
# true k
k <- 3
n <- nrow(X)
p <- ncol(X)
A <- matrix(as.numeric(!R), nrow = n, ncol = p)
ps <- rowSums(A)
# Unique patterns of missingness.
R.unique <- unique(R, MARGIN = 1)
# Missingness pattern of each observation.
miss.grp <- apply(R, 1, MixtClust:::row_match, R.unique)
# Selection indices for each pattern of missingness.
Ru <- 1*!R.unique
# Initial values
init <- MixtClust:::get.init.val(Y, R, k, FALSE, FALSE)
mus <- init$mu; Sigmas <- init$Sigma; pis <- init$pi; nus <- init$nu


###########################################################################
context("Density")

test_that("Mahalanobis distance", {
  expect_equivalent(stats::mahalanobis(x = Y, center = mus[k,], cov = Sigmas[,,k]), 
                    as.numeric(mahalanobis(x = Y, mu = mus[k,], sigma = Sigmas[,,k])))
})

test_that("Array-cube conversion", {
  expect_equivalent(Sigmas, MixtClust:::to_array(as.numeric(Sigmas), k, p))
})

test_that("mvt t density", {
  expect_equivalent(as.numeric(dMVT(Y, mus[k,], Sigmas[,,k], nus[k])), 
                    mvnfast::dmvt(Y, mus[k,], Sigmas[,,k], nus[k]))
})

test_that("Observed mvt t density", {
  expect_equivalent(as.numeric(h(Y, mus[k,], Sigmas[,,k], nus[k], miss.grp, Ru)),
                    sapply(1:n, function(i) mvnfast::dmvt(Y[i,!R[i,]], mus[k,!R[i,]], Sigmas[!R[i,],!R[i,],k], nus[k])))
})


###########################################################################
context("E steps")

ZZ <- MixtClust:::up_Z(Y, mus, Sigmas, nus, pis, miss.grp, Ru)
test_that("E-step Z", {
  expect_equivalent(rowSums(ZZ), 
                    rep(1, n))
})

WW <- MixtClust:::up_W(Y, mus, Sigmas, nus, miss.grp, Ru)


###########################################################################
context("M steps")

test_that("M-step mixing proportions", {
  expect_equal(sum(MixtClust:::up_pi(ZZ)), 1)
})

test_that("M-step mu (marginalization)", {
  expect_equivalent(MixtClust:::up_mu(Y, ZZ, WW, A),
                    MixtClust:::up_muR(Y, ZZ, WW, A))
})

M <- length(unique(miss.grp))
SOiOEOO <- lapply(1:k, function(kk) MixtClust:::SOiOEOOk(Sigmas[,,kk], Ru))
xhat <- lapply(1:k, function(kk) MixtClust:::xhatk(Y, mus[kk,], miss.grp, M, SOiOEOO[[kk]]))
mus3 <- MixtClust:::up_mu_Lin(p, ZZ, WW, xhat)
test_that("Lin's method imputation", {
  expect_identical(sapply(1:k, function(kk) all.equal(Y[!R], xhat[[kk]][!R])), 
                   rep(TRUE, k))
})

test_that("M-step sigma (marginalization)", {
  expect_equivalent(MixtClust:::up_Sigma(Y, ZZ, WW, mus, A, FALSE),
                    MixtClust:::up_SigmaR(Y, ZZ, WW, mus, A, FALSE))
  
})

SigmasC <- MixtClust:::up_Sigma(Y, ZZ, WW, mus, A, TRUE)
test_that("M-step sigma constrained (marginalization)", {
  expect_identical(sapply(2:k, function(kk) all.equal(SigmasC[,,1], SigmasC[,,kk])), 
                   rep(TRUE, k-1))
})


test_that("Ruegree of freedom constraints", {
  expect_true(length(unique(MixtClust:::up_nu(ZZ, WW, nus, ps, TRUE, FALSE))) == 1)
  expect_true(length(unique(MixtClust:::up_nu(ZZ, WW, nus, ps, TRUE, TRUE))) == 1)
  
})


###########################################################################
context("Agreement on complete cases")
# Compare marginalization to full EM updates for mu, sigma w/ CC (full data d)
gd <- rep(1, n); Rud <- matrix(1, 1, p); Ad <- matrix(1, n, p)
ZZd <- MixtClust:::up_Z(as.matrix(d), mus, Sigmas, nus, pis, gd, Rud)
WWd <- MixtClust:::up_W(as.matrix(d), mus, Sigmas, nus, gd, Rud)
mud1 <- MixtClust:::up_mu(as.matrix(d), ZZd, WWd, Ad)
SOiOEOOd <- lapply(1:k, function(kk) MixtClust:::SOiOEOOk(Sigmas[,,kk], Rud))
xhatd <- lapply(1:k, function(kk) MixtClust:::xhatk(as.matrix(d), mus[kk,], gd, 1, SOiOEOOd[[kk]]))
mud2 <- MixtClust:::up_mu_Lin(p, ZZd, WWd, xhatd)
Sigmasd1 <- MixtClust:::up_Sigma(as.matrix(d), ZZd, WWd, mud1, Ad, FALSE)
Sigmasd2 <- sapply(1:k, 
                   function(kk) MixtClust:::up_Sigmak_Lin(1, ZZd[,kk], WWd[,kk], mud2[kk,], Sigmas[,,kk], xhatd[[kk]], gd, SOiOEOOd[[kk]]), 
                   simplify = 'array')
test_that("Mu updates on CC: fullEM vs marginalization", {
  expect_equivalent(mud1, mud2)
})
test_that("Sigma updates on CC: fullEM vs marginalization", {
  expect_equivalent(Sigmasd1, Sigmasd2)
})
seed <- round(runif(1)*1e3, 0)
# Complete Case (marginalization code)
set.seed(seed)
ans.marg.CC <- MixtClust(X[rowSums(R)==0,], nclusters = k, emEM.args = list(nstarts=k*10, em.iter=5, nbest=1))
lln <- ans.marg.CC$loglikn + MixtClust:::loglikelihood(X[rowSums(R)>0,], ans.marg.CC$estimates)
lln2 <- MixtClust:::loglikelihood(X, ans.marg.CC$estimates)
# Complete Case (deletion)
set.seed(seed)
ans.del <- MixtClust(X, nclusters = k, emEM.args = list(nstarts=k*10, em.iter=5, nbest=1), method = 'deletion')
# Complete Case (full EM code)
set.seed(seed)
ans.fullEM.CC <- MixtClust(X[rowSums(R)==0,], nclusters = k, emEM.args = list(nstarts=k*10, em.iter=5, nbest=1), method = 'fullEM')
test_that("Main function on CC: fullEM vs marginalization", {
  expect_equivalent(ans.marg.CC$estimates, ans.fullEM.CC$estimates)
  expect_equivalent(ans.marg.CC$class, ans.fullEM.CC$class)
})
test_that("Main function: marginalization on CC vs deletion", {
  expect_equivalent(ans.marg.CC$estimates, ans.del$estimates)
  expect_equivalent(ans.marg.CC$class, ans.del$class[rowSums(R)==0])
})
test_that("Loglik calculation", {
  expect_equivalent(ans.del$loglikn, lln)
  expect_equivalent(lln2, lln)
})