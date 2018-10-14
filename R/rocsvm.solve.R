#' Finding Lagrangian multipliers of ROC-SVM by Qudratic Programming
#'
#' Computes the Lagrangian multipliers(\code{alpha}), which are solutions of ROC-SVM using Quadratic Programming.
#'
#' @param K The kernelized matrix, i.e., K< .,. >.
#' @param lambda The regularization parameter that users want in ROC-SVM model.
#' @param rho A positive constant (default : 1)
#' @param eps Adjustment computing errors (default : 1e-08)
#' @author Seung Jun Shin, Do Hyun Kim
#' @seealso \code{\link{rocsvm.path}}
#' @examples
#' \donttest{
#' n <- 30
#' p <- 2
#' delta <- 1
#' set.seed(309)
#' y <- c(rep(1, n/2), rep(-1, n/2))
#' x <- matrix(0, n, p)
#' for (i in 1:n){
#'  if (y[i] == 1) {
#'  x[i,] <- rnorm(p, -delta, 1)
#'  } else {
#'  x[i,] <- rnorm(p, delta, 1)
#'   }
#'  }
#'
#' K <- radial.kernel(x,x)
#' rocsvm.solve(K, lambda = 1, rho = 1) }
#' @export
#' @importFrom quadprog solve.QP
rocsvm.solve <- function(K, lambda, rho = 1, eps = 1.0e-8)
{
  N <- dim(K)[1]
  obj <- solve.QP(Dmat = (K + diag(rep(eps, N)))/lambda,
                  dvec = rep(rho, N),
                  Amat = cbind(diag(N), -diag(N)),
                  bvec = c(rep(0, N), rep(-1, N)))
  value <- list(alpha = obj$solution, value = obj$value)
  return(value)
}

