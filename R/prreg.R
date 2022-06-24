#'Hybrid predictive recursion-EM algorithm for regression
#'
#'Fits a robust linear regression model
#'
#'@param X Matrix of covariate values
#'@param Y Vector of response values
#'@param perm A matrix of user-defined permutation sequence containing \code{nperm} rows and number of columns equal to number of data points in X
#'@param get.ci Logical. Should confidence intervals be calculated for the coefficients? Default is [FALSE]
#'@return A list with the following elements,
#' \itemize{
#'   \item b - Estimated regression coefficients
#'   \item Yhat - Predicted responses at X
#'   \item ci - Confidence intervals for b
#'   \item f - Permutation-averaged mixing density estimate
#'   \item UU - Vector of support for f
#'   \item wt - Final vector of weights in the regression fit
#'   \item lik - Permutation-averaged log-likelihood function L
#' }
#'
#'@examples
#'# Belgian phone call data in Example 1 of the paper
#'phone.example <- function() {
#'require(MASS)
#'data(phones)
#'Y <- phones$calls
#'X <- phones$year + 1900
#'o <- prem.reg(X=cbind(1 + 0 * Y, X), Y=Y)
#'plot(x=X, y=Y, xlab="Year", ylab="Calls (in millions)")
#'abline(a=o$b[1], b=o$b[2])
#'abline(lm(Y ~ X), lty=2)
#'legend(x="topleft", inset=0.05, lty=1:2, c("PREM", "LS"))
#'return(list(b=o$b, f=o$f, UU=o$UU, w=o$wt, lik=o$lik))
#'}
#'
#'#TO RUN
#'phone.example()
#'@source A first version of the code is available at \url{https://www4.stat.ncsu.edu/~rmartin/Codes/prreg.R}
#'@references
#'Martin, Ryan, and Zhen Han. "A semiparametric scale-mixture regression model
#'and predictive recursion maximum likelihood."
#'Computational Statistics & Data Analysis 94 (2016): 75-85.
#'
#'@seealso [pr()]
#'
#'@export
prem.reg <- function(X, Y, perm, get.ci=FALSE) {

  tol <- 1e-05
  ols <- lm(Y ~ X - 1)
  sigma <- summary(ols)$sigma
  UUmax <- max(5 * sigma, 50)
  UU <- seq(1e-05, UUmax, len=201)
  K <- function(z, u) dnorm(z, 0, u)
  if(missing(perm)) {

    n <- length(Y)
    nperm <- 25
    perm <- matrix(0, ncol=nperm, nrow=n)
    perm[,1] <- 1:n
    for(j in 2:nperm) perm[,j] <- sample(n)

  }
  else nperm <- ncol(perm)
  b.old <- as.numeric(ols$coefficients)
  f <- 1 + 0 * UU
  f <- f / int(f, UU)
  lik <- c()
  repeat {

    Z <- as.numeric(Y - X %*% b.old)
    o <- pr.reg(X=Z, d=K, f0=f, U=UU, perm=perm)
    lik <- c(lik, -o$L)
    ww <- o$wt
    b <- as.numeric(lm(Y ~ X - 1, weights=ww)$coefficients)
    if(sum(abs(b - b.old)) < tol) break else b.old <- b

  }
  Z <- as.numeric(Y - X %*% b)
  o <- pr.reg(X=Z, d=K, f0=f, U=UU, perm=perm)
  lik <- c(lik, -o$L)
  ww <- o$wt
  if(get.ci) {

    if(!exists("hessian")) library(numDeriv)   #requires the R package "numDeriv"
    prml <- function(bb) {

      Z <- as.numeric(Y - X %*% bb)
      L <- pr.reg(X=Z, d=K, f0=f, U=UU, perm=perm)$L
      return(L)

    }
    J <- solve(hessian(prml, b))
    ci <- b + 1.96 * outer(sqrt(diag(J)), c(-1, 1))

  } else ci <- NULL
  return(list(b=b, Yhat=as.numeric(X %*% b), ci=ci, f=o$f, UU=UU, wt=ww, lik=lik))

}
#'
#'@export
pr.reg <- function(X, d, U, f0, w, perm) {

  n <- length(X)
  if(missing(f0)) { f0 <- 1 + 0 * U; f0 <- f0 / int(f0, U) }
  if(missing(w)) w <- function(i) 1 / (i + 1)
  N <- ncol(perm)
  f.avg <- 0 * U
  L.avg <- 0
  wt <- 0 * X
  for(j in 1:N) {

    f <- f0
    L <- 0
    for(i in 1:n) {

      num <- d(X[perm[i, j]], U) * f
      den <- int(num, U)
      wt[perm[i, j]] <- wt[perm[i, j]] + int(num / U**2, U) / den
      L <- L + log(den)
      f <- (1 - w(i)) * f + w(i) * num / den

    }
    f.avg <- f.avg + f / N
    L.avg <- L.avg + L / N

  }
  return(list(f=f.avg, L=-L.avg, wt=wt / N))

}


