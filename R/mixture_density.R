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
#' \item A matrix where each row is a random sample from f0 for multivariate f
#' }
#'@param d Parametric kernel function of form d(x,U)
#'@param Xsup Mixture density support, where Xsup is a
#' \itemize{
#' \item Equispaced and odd length vector for a univariate mixture
#' \item Matrix of two columns, where each column is the support for the corresponding variate for a bivariate mixture
#' \item A matrix where each row is a support point for a multivariate mixture
#' }
#'@return Vector containing mixture density estimate values at Xsup
#'
#'@seealso [pr()]

# Usage
#'@export
mixture_density <- function(f, U, d, Xsup,...){
  U = as.matrix(U)
  Xsup = as.matrix(Xsup)
  if(ncol(Xsup)==2) {
    Xsup = as.matrix(expand.grid(Xsup[,1],Xsup[,2]))
  }
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
  mix = list(Xsup = Xsup, mix = mix)
  class(mix) = "mixture"
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
  num = sapply(U, d, x = x)*f
  return(int(num, U))
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
  num = d(x=x, u=U.l)*f
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
#'@export
single.mix.multi = function(x, D, U, d, t,...){
  ans = apply(U, 1, d, x = x)
  return((1/t)*sum(ans*D))
}

#'Plots the PR mixture density function
#'@export
plot.mixture = function(obj){
Xsup = obj$Xsup
dx = ncol(Xsup)
if(dx==1){
  plot(Xsup, obj$mix, xlab = "X", ylab = "m", type = "l", main = "Estimated mixture density")
  }
else if(dx==2) {
  s1 = unique(Xsup[,1]) ; s2=unique(Xsup[,2])
  m.matrix = matrix(obj$mix, length(s1), length(s2), byrow = TRUE)
  ContourFunctions::gcf_grid(s2, s1, m.matrix, mainminmax = FALSE, color.palette = function(x) rev(gray((1:x)/x)), bar = TRUE)
}
else {
  print("Plot yet to be decided")
}
}
