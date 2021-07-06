#'
#' @name MixtClust
#'
#' @title Robust Clustering for Complete and Incomplete Data
#' 
#' @description Robust clustering, including handling of incomplete data, using the EM algorithm for 
#'   finite mixtures of multivariate t distributions
#'
#' @details Model-based clustering using finite mixtures of t distributions, with handling of incomplete data using either marginalization or the EM algorithm.
#'   
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @importFrom stats runif
#' @importFrom stats cov.wt
#' @importFrom stats kmeans
#' @importFrom kmmeans kmmeans
#' @importFrom utils tail
#' 
#' @useDynLib MixtClust, .registration = TRUE
#' 
NULL