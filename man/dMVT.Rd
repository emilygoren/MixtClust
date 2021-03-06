% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{dMVT}
\alias{dMVT}
\title{Multivariate t distribution}
\usage{
dMVT(x, mu, sigma, nu, logans = FALSE, ischol = FALSE)
}
\arguments{
\item{x}{Numeric. A vector (of length \eqn{p}) or matrix (with \eqn{p}
columns).}

\item{mu}{Numeric. A vector of length \eqn{p} representing the mean.}

\item{sigma}{Numeric. A \eqn{p \times p} non-negative definite matrix (or
its Cholesky decomposition).}

\item{nu}{Numeric. A positive scalar representing the degrees of freedom.}

\item{logans}{Logical. If \code{TRUE}, the log density is returned.}

\item{ischol}{Logical. Set to \code{TRUE} if \code{sigma} is provided as a
Cholesky decomposition.}
}
\value{
The multivariate t density for all rows in \code{x} using degrees of
  freedom \code{nu}, mean vector \code{mu}, and covariance matrix
  \code{sigma}.
}
\description{
Compute the density of the multivariate t distribution (fast
  implementation).
}
\author{
Emily Goren, \email{emily.goren@gmail.com}
}
