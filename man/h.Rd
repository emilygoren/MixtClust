% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{h}
\alias{h}
\title{Multivariate t distribution with missing data}
\usage{
h(x, mu, sigma, nu, grp, Ru)
}
\arguments{
\item{x}{Numeric. A vector (of length \eqn{p}) or matrix (with \eqn{p}
columns). Missing values for row \eqn{i} correspondind to unique pattern
of missingness \eqn{m} are indicated by the \eqn{m^{th}} element of the
selection index list \code{D}.}

\item{mu}{Numeric. A vector of length \eqn{p} representing the mean.}

\item{sigma}{Numeric. A \eqn{p \times p} non-negative definite matrix (or
its Cholesky decomposition).}

\item{nu}{Numeric. A postitive scalar.}

\item{grp}{Numeric. A vector of length \eqn{n} with elements indicating the
which row of \code{Ru} to use as the observed coordinates for each of
the \eqn{n} observations.}

\item{Ru}{Binary matrix. Each row corresponds to a unique pattern of missingness,
where \eqn{1} indicates observed and \eqn{0} indicates missing coordinates.}
}
\value{
The multivariate t density for all rows in \code{x} using degrees of
  freedom \code{nu}, mean vector \code{mu}, covariance matrix \code{sigma},
  and observed coordinates specified by rows in \code{Ru}.
}
\description{
Compute the marginal density of the multivariate t distribution
  for the observed coordinates.
}
\author{
Emily Goren, \email{emily.goren@gmail.com}
}
