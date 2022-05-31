# Description
#' Mixture density estimate 
#'
#' Calculates the mixture density at a user-defined grid of points given a mixing density
#'
#'@param f Vector of mixing density values or vector of particle approximation over support U
#'@param U Mixing density support, where U is a
#' \itemize{
#' \item Equispaced and odd length vector for univariate f 
#' \item Matrix of two columns, where each column is the support for the corresponding variate for bivariate f
#' \item A matrix where each row is a random sample from f0
#' }          
#'@param d Parametric kernel function of form d(x,u)
#'@param Xsup Vector or matrix (each row representing one support point) at which mixture density is calculated 
#'
#'@return Vector containing mixture density estimate values at Xsup
#'
#'@examples
#'@seealso [pr()]

# Usage
mixture_density <- function(f, U, d, Xsup,...){
  U = as.matrix(U)
  t <- nrow(U)
  du <- ncol(U)
  if(du==1){
    mix = sapply(Xsup, single.mix1, f = f, U = U, d = d,...)
  } 
  else if(du==2) {
    mix = apply(Xsup, 1, single.mix2, f = f, U = U, d = d, t = t,...)
  } else {
    mix = apply(Xsup, 1, single.mix.multi, U = U, D = f, t = t, d = d,...)
  }
  return(mix)
}

#' Mixture density at point x for unidimensional U
#' 
#' @param x Single value at which the mixture density is to be calculated 
#' @param f See documentation for mixture_density
#' @param U See documentation for mixture_density
#' @param d See documentation for mixture_density
#' 
#' @return Mixture density calculated at x
#'
#'@seealso [pr()], [mixture_density()]
#'@export
single.mix1 = function(x, f, U, d,...){
  return(int(d(x, U,...)*f, U))
}

#' Mixture density at point x for bidimensional U
#' 
#' @param x Single value at which the mixture density is to be calculated
#' @param f See documentation for mixture_density
#' @param U See documentation for mixture_density
#' @param d See documentation for mixture_density
#' @param t Number of rows of U 
#' 
#' @return Mixture density calculated at x
#'@seealso [pr()], [mixture_density()]
#'@export
single.mix2 = function(x, f, U, d, t,...){
  U.l = as.matrix(expand.grid(U[,1],U[,2]))
  num = d(x, U.l,...)*f
  f_matrix = matrix(num, nrow = t, ncol = t, byrow=TRUE)
  return(simp.int2(U[,1], U[,2], f_matrix))
}

#' Mixture density at point x for multidimensional U
#' 
#' @param x Single value at which the mixture density is to be calculated
#' @param D Mixing density weights at U
#' @param U See documentation for mixture_density
#' @param d See documentation for mixture_density
#' @param t Number of rows of U 
#' 
#' @return Mixture density calculated at x
#'
#'@seealso [pr()], [mixture_density()]
single.mix.multi = function(x, D, U, d, t,...){
  return((1/t)*sum(d(x, U,...)*D))
}
