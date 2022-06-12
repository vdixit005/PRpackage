#'Monotone density estimation using PR
#'
#'Estimates a monotone density via the PR algorithm
#'
#'@param X Vector of data values
#'@param nperm Number of permutations for the PR estimate
#'
#' @return A list with the following elements,
#' \itemize{
#'   \item m - Estimated monotone density
#'   \item Xsup - Corresponding support
#' }
#'
#'@examples
#'# To estimate a monotone density using PR
#'# Generate data
#'n = 200
#'X = abs(rnorm(n))
#'# Estimate
#'m = monotone(X)
#'# Plot with the truth
#'truth = 2*dnorm(m$Xsup)
#'plot(m)
#'lines(m$Xsup, m$m, col = "red")

#' @references
#' Dixit, Vaidehi, and Ryan Martin. "Revisiting consistency of a
#' recursive estimator of mixing distributions."
#' arXiv preprint arXiv:2110.02465 (2021).
#'@export
monotone = function(X, nperm = 25){
n = length(X)
U = seq(10^-5, max(X), length.out = 101)
f0 = dunif(U, min = min(U), max = max(U))
Xsup = seq(10^-20, max(X), length.out = 101)
f = pr(X = X, d = d.mon, U = U, f0 = f0, nperm = nperm)
m = mixture_density(f = f$f, U = U, d = d.mon, Xsup = Xsup)
return(m)
}
#'
d.mon = function(x,u){
  dunif(x, min = 0, max = u)
}

