## usethis namespace: start
#' @useDynLib PredictiveRecursion, .registration = TRUE
## usethis namespace: end
NULL
## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
#'
#' @title Predictive recursion (PR) algorithm
#'
#' @description Executes the PR algorithm
#'
#' @param X Vector or matrix of data
#' @param d Parametric kernel function of form d(x,u)
#' @param U Mixing density support, where U is a
#' \itemize{
#' \item Equispaced and odd length vector as the support for univariate f
#' \item Matrix of two columns, where each column is the support for the corresponding variate for bivariate f
#' \item A matrix where each row is a random sample from f0 for multivariate f
#' }
#'See 'Details'
#' @param f0 A vector of initial guess for f
#' @param w A weight function defining the weight sequence
#' @param nperm Number of permutations to be used
#' @param perm A matrix of user-defined permutation sequence containing \code{nperm} rows and number of columns equal to number of data points in X
#' @param ... Additional arguments to d
#'
#' @return A list with the following elements,
#' \itemize{
#'   \item f - Permutation-averaged mixing density estimate
#'   \item L - Permutation-averaged negative log-likelihood function L
#'   \item D - Permutation-averaged weights D at U (for multivariate case only)
#' }
#' @section Details:
#' The PR algorithm is used to estimate mixing density \eqn{f} given data \eqn{X_1, \ldots, X_n} from mixture density \eqn{m} with a known kernel \eqn{k}, where
#' \deqn{m(x) = \int k(x \mid u) f(u) du}
#'
#' This is a sequential algorithm which starts with a user-defined guess \eqn{f0} and returns the estimate \eqn{fn} in \eqn{n} steps.
#' If \eqn{f(u)} is univariate then this integration is approximated by numerical integration; where the user specifies an equispaced and odd length vector \eqn{U} as the support.
#' If \eqn{f(u)} is bivariate i.e. \eqn{u = (u_1, u_2)} then this integration is approximated with a quadrature rule. For this a two-dimensional grid is created by taking the first column of \eqn{U} as support for \eqn{u_1} and second column of \eqn{U} as support for \eqn{u_2}.
#' Lastly, if \eqn{f(u)} is multivariate i.e. \eqn{u = (u_1,...,u_d)}, then the integration is approximated as a weighted sum over random observations \eqn{U}. In this case, each row of \eqn{U} is a vector of length \eqn{d} and is a random observation from \eqn{f0}.
#' @examples
#' # Univariate implementation
#'U = seq(0, 1, length.out = 101)
#'f0 = dunif(U, min = 0, max = 1)
#'
#'# Generate data from a univariate mixture density
#'u = rbeta(n = 1000, 2, 10)
#'x = rnorm(n = 1000, mean = u, sd = 0.1)
#'
#'# PR estimate
#'ans1 = pr(X = x, d = dnorm, U = U, f0 = f0, sd = 0.1)
#'plot(ans1)
#'lines(U, dbeta(U, 2,10), col = "red")
#'
#'# Mixture density estimate
#'Xsup = seq(-0.3, 1.3, length.out = 101)
#'m1 = mixture_density(f = ans1$f, U = U, d = dnorm, Xsup = Xsup, sd = 0.1)
#'plot(m1)
#'############################################################################################################
#'# Bivariate implementation
#'U = cbind(seq(0, 1, length.out = 101), seq(10^-5, 1, length.out = 101))
#'f0 = rep(1, 101*101)
#'d2 = function(x,u){
#'  return(dnorm(x, mean = u[,1], sd = u[,2]))
#'}
#'
#'# Generate data from a bivariate mixture density
#'u = cbind(rbeta(n = 2000, 2, 10), rbeta(n = 2000, 4, 4))
#'x = rnorm(n = 2000, mean = u[,1], sd = u[,2])
#'
#'# PR estimate
#'ans2 = pr(X = x, d = d2, U = U, f0 = f0)
#'plot(ans2)
#'
#U.l = as.matrix(expand.grid(U[,1],U[,2]))
#f.truth = dbeta(U.l[, 1], 2, 10) * dbeta(U.l[, 2], 4, 4)
#ContourFunctions::gcf_grid(U[,2], U[,1], matrix(f.truth, 101, 101, byrow = TRUE))
#'# Mixture density estimate
#'Xsup = as.matrix(seq(-1.5, 2.8, length.out = 101))
#'m2 = mixture_density(f = ans2$f, U = U, d = d2, Xsup = Xsup)
#'plot(m2)
#m.truth = mixture_density(f.truth, U, d = d2, Xsup)
#'
#'##########################################################################################
#'
#'# Multivariate implementation
#'require(sn)
#'d3 = function(x,u){
#'ans = dsn(x, xi = u[1], omega = u[2], alpha = u[3])
#'return(ans)
#'}
#'
#'# Generate particles from initial guess f0
#'t = 3000
#'U = cbind(runif(t, min = 0, max = 1), runif(t, min = 0, max = 1), runif(t, min = -4, max = 4))
#'
#'# PR estimate
#'ans3 = pr(X = x, d = d3, U = U)
#'plot(ans3)
#'
#'# Mixture density estimate
#'Xsup = as.matrix(seq(-1.5, 2.8, length.out = 101))
#'m3 = mixture_density(f = ans3$D, U = U, d = d3, Xsup = Xsup)
#'plot(m3)
#'
#'@source A first version of the code is available at \url{https://www4.stat.ncsu.edu/~rmartin/Codes/pr.R}
#'@references
#'Dixit, Vaidehi, and Ryan Martin. "Estimating a mixing
#'distribution on the sphere using predictive recursion."
#'Sankhya B (2022): 1-31.
#'
#'Dixit, Vaidehi, and Ryan Martin. "A PRticle filter algorithm for
#' nonparametric estimation of multivariate mixing distributions."
#' arXiv preprint arXiv:2204.01646 (2022).
#'
#'############################################################################################################################
#'@export
pr <- function(X, d, U, f0, w, nperm = 1, perm = NULL,...) {
  X = as.matrix(X)
  U = as.matrix(U)
  n <- nrow(X)
  t <- nrow(U)
  du <- ncol(U)
  if(missing(f0)) f0 <- 1 + 0 * U[,1]
  if(missing(w)) w <- function(i) 1 / (i + 1)^{0.67}
  N <- nperm
  if(N == 1 && is.null(perm)) {
    perm <- matrix(0, n, N)
    perm[,1] <- 1:n}
  else if(N > 1 && is.null(perm)) {
    perm <- matrix(0, n, N)
    perm[,1] <- 1:n
    for(j in 2:N) perm[,j] <- sample(n)

  }
  f.avg = 0 * f0
  L.avg = 0
  if(du==1){
    f0 = f0 / int(f0, U)
    for(j in 1:N) {
      f = f0
      L = 0
      x <- X[perm[,j],]
      for(i in 1:n) {
        num <- d(x[i], U,...) * f
        den <- int(num, U)
        L <- L + log(den)
        f <- (1 - w(i)) * f + w(i) * num / den
      }
      f.avg <- (j - 1) * f.avg / j + f / j
      L.avg <- (j - 1) * L.avg / j + L / j
    }
    ans = list(U = U, f=f.avg, L=-L.avg)
    class(ans) = "pr"
    return(ans)
  }
  else if(du==2){
    f0 = f0 / simp.int2(U[,2], U[,1], matrix(f0, t, t, byrow=TRUE))
    U.l = as.matrix(expand.grid(U[,1],U[,2]))
    for(j in 1:N) {
      f = f0
      L = 0
      x <- as.matrix(X[perm[,j], ])
      for(i in 1:n){
        num = d(x[i,], U.l,...) * f
        f_matrix = matrix(num, nrow = t, ncol = t, byrow=TRUE)
        den = simp.int2(U[,2], U[,1], f_matrix)
        L <- L + log(den)
        f <- (1 - w(i)) * f + w(i) * num / den
      }
      f.avg <- (j - 1) * f.avg / j + f / j
      L.avg <- (j - 1) * L.avg / j + L / j
    }
    ans = list(U = U, f=f.avg, L=-L.avg)
    class(ans) = "pr"
    return(ans)
  }
  else {
    # D.avg = 0 * f0
    # for(j in 1:N) {
    #   f = f0
    #   L = 0
    #   x <- as.matrix(X[perm[,j], ])
    #   D = 1
    #   for(i in 1:n) {
    #     num <- d(x[i,], U,...) * f
    #     den <- (1/t)*sum(d(x[i,], U,...)*D)
    #     L <- L + log(den)
    #     f <- (1 - w(i)) * f + w(i) * num / den
    #     D = D*(1 + w(i)*(d(x[i,], U,...)/den - 1))
    #
    #   }
    #   # Add D.avg
    #   f.avg <- (j - 1) * f.avg / j + f / j
    #   L.avg <- (j - 1) * L.avg / j + L / j
    #   D.avg <- (j - 1) * D.avg / j + D / j
    w.vec = w(1:n)
    cpp.final = pr_cpp(f0, U, X, d, N, w.vec)
    ans = list(U = U, f=cpp.final$f, L=-cpp.final$L, D = cpp.final$D)
    class(ans) = "pr"
    return(ans)
  }
}

#'#Plots the PR density function
#'@export
plot.pr = function(obj){
  U = as.matrix(obj$U)
  du = ncol(U)
  if(du==1){
  plot(U, obj$f, xlab = "U", ylab = "f", type = "l", main = "Estimated mixing density")
  }
  else if(du==2) {
  f.matrix = matrix(obj$f, nrow(U), nrow(U), byrow = TRUE)
  ContourFunctions::gcf_grid(U[,2], U[,1], f.matrix, mainminmax = FALSE, color.palette = function(x) rev(gray((1:x)/x)), bar = TRUE)
  }
  else {
    print("Plot yet to be decided")
  }
}

########## Adding a C++ function to replace the multivariate sum in pr.R



