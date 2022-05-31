#' Auxiliary univariate numerical integration
#' 
#' @param f Vector of values of function f at x
#' @param x Equispaced and odd length vector
#' @param tol Error tolerance
#' 
#' @return Scalar value of numerical integration of function f over x
int <- function(f, x, tol=1e-10) {
  
  n <- length(x)
  if(n == 1) return(0)
  simp <- (sum(range(x[-1] - x[-n]) * c(-1, 1)) < tol)
  if(!simp) out <- sum((f[-1] + f[-n]) * (x[-1] - x[-n])) / 2
  else out <- ((x[2] - x[1]) / 3) * (sum(f) + sum(f[2:(n-1)]) + 2 * sum(f[seq(2, n-1, 2)]) )
  return(out)
  
}

#' Auxiliary bivariate numerical integration
#' 
#' @param x Equispaced vector
#' @param y Equispaced vector
#' @param Fxy Matrix of dimension dim(x) * dim(y) that contains values of f(x,y)
#' @param tol Error tolerance
#' 
#' @return Scalar value of numerical integration of function f over x
simp.int2 <- function(x, y, Fxy, tol=1e-10) {
  
  n <- length(x); if(n %% 2 == 0) stop("need an odd number of grid points")
  m <- length(y); if(m %% 2 == 0) stop("need an odd number of grid points")
  
  check <- (sum(range(x[-1]-x[-n])*c(-1,1)) < tol) && (sum(range(y[-1]-y[-m])*c(-1,1)) < tol)
  if(!check) stop("xy-grid must be equi-spaced")
  
  h <- x[2] - x[1]  # spacings for x (assumed to be equi-spaced)
  k <- y[2] - y[1]  # spacings for y
  
  u <- c(1, rep(c(4,2), (n-3)/2), 4, 1)
  v <- c(1, rep(c(4,2), (m-3)/2), 4, 1)
  M <- outer(u, v, '*')
  
  return(h * k * sum(M * Fxy) / 9)
  
}