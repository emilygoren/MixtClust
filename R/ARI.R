

#' Adjusted Rand Index
#'
#' @param class1 A \eqn{n}-vector of class labels.
#' @param class2 A \eqn{n}-vector of class lables.
#' 
#' @details Compute the adjusted Rand Index comparing the two 
#' classifications, \code{class1} and \code{class2}.
#'
#' @references 
#' Lawrence Hubert & Phipps Arabie, 1985. "Comparing partitions," Journal of
#' Classification, 2(1): 193-218.
#'
#' @return A scalar value, taking the value 1 if the two classifications are
#' in perfect agreement. Under random classification, the adjusted Rand Index is 0.
#' 
#' @author Emily Goren, \email{emily.goren@gmail.com}
#'
#' @export
adjRand <- function(class1, class2) {
  x <- as.integer(class1)
  y <- as.integer(class2)
  if (all(x == y))
    return(1)
  n <- length(x)
  if (n != length(y)) 
    stop("Please supply arguments with the same length.")
  contab <- table(x, y)
  index <- sum(choose(contab, 2))
  a <- rowSums(contab)
  b <- colSums(contab)
  A <- sum(choose(a, 2)) 
  B <- sum(choose(b, 2))
  Eindex <- A * B / choose(n, 2)
  maxindex <- (A + B) / 2
  out <- (index - Eindex) / (maxindex - Eindex)
  if (is.nan(out)) out <- 0
  return(out)
}


