

#' Classify New Observations
#'
#' @description Classify new, possibly incomplete, observations arising from a finite mixture of multivariate
#' t distributions.
#' 
#' @param newobs A matrix with new observations (rows), \eqn{p} columns
#'   (dimensions), and missing entries set to \code{NA}.
#' @param params A list of parmaters defining a finite mixture of
#' multivariate t distributions, 
#' usually obtained from \code{\link{MixtClust}} (see details).
#' 
#' @details Classify new observations according to the finite mixture of t distirbutions
#' specified by the parameter values in \eqn{params}, a named list with elements:
#' \itemize{ 
#'   \item{"pi"}{Mixing proportions. A vector of length \eqn{K} that
#' sums to one.} 
#'   \item{"nu"}{Degrees of freedom. A vector of length \eqn{K} with
#' entries at least equal to three (thus requiring the existance of the first
#' two moments.)} 
#'   \item{"mu"}{Locations. A \eqn{K \times p} matrix, where the
#' \eqn{k}-th row is the location \eqn{\mu_k \in R^p} for cluster
#' \eqn{k}.} 
#'   \item{"Sigma"}{Dispersions. A \eqn{p \times p \times K} array, where
#' the \eqn{k}-th slice is the \eqn{p \times p} positive-definite disperion
#' matrix \eqn{\Sigma_k} for cluster \eqn{k}.} }
#'
#' @return A vector classify each observation to the cluster \eqn{1, \dots, K}
#' with the highest posterior probability.
#' 
#' @author Emily Goren, \email{emily.goren@gmail.com} based on modifications of
#' code by Ranjan Maitra.
#'
#' @export
classify <- function(newobs, params) {
  if (is.null(dim(newobs))) {
    x <- t(as.matrix(newobs))
  } else {
    x <- as.matrix(newobs)
  }
  prior <- params$pi
  df <- params$nu
  K <- length(prior)
  if (ncol(params$mu) != ncol(x))
    stop("Dimension of newobs and mu do not match")
  if (any(dim(params$Sigma)[1:2] != ncol(x)))
    stop("Dimension of newobs and Sigma do not match")
  out <- apply(x, 1, function(xi) {
    rmv <- !is.na(xi)
    if (sum(!rmv) == length(xi)) {
      ans <- NA
    } else {
      x.tr <- t(matrix(xi[rmv]))
      mean.tr <- params$mu[, rmv]
      sigma.tr <- params$Sigma[rmv, rmv, ]
      if (sum(rmv) == 1) {
        postr <- sapply(1:K, function(k) prior[k] * dMVT(x.tr, as.matrix(mean.tr[k]), as.matrix(sigma.tr[k]), df[k]))
      } else {
        postr <- sapply(1:K, function(k) prior[k] * dMVT(x.tr, mean.tr[k,], sigma.tr[,,k], df[k]))
      }
      postr <- postr/sum(postr)
      ans <- which.max(postr)
    }
    return(ans)
  })
  return(t(out))
}