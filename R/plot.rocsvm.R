#' Plot the rocsvm.path, solution paths of ROC-SVM as a function of lambda
#'
#' produces a plot of the ROC-SVM \code{lambda} path.
#'
#' @param x The rocsvm path object
#' @param ... Generic compatibility
#' @return The entire solution path of ROC-SVM solution as a function of \code{lambda}.
#' @author Seung Jun Shin, Do Hyun Kim
#' @seealso \code{\link{rocsvm.path}}
#' @examples
#' # The 'obj' comes from an example description of rocsvm.path()
#' library(rocsvm.path)
#'
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
#' rho = 1
#' kernel = radial.kernel
#' param.kernel  = 1/ncol(x)
#' prop = 0.1
#'
#' obj <- rocsvm.path(x, y, rho, kernel, param.kernel, prop)
#' plot(obj)
#' # or plot.rocsvm(obj, lty = 2, lwd = 2, col = 2)
#'
#' @export
#' @export plot.rocsvm
#' @importFrom graphics abline lines par plot
#'
plot.rocsvm <- function(x, ...){
  lambda <- x$lambda
  alpha  <- x$alpha
  N <- nrow(alpha)

  par(mfrow = c(1,1), mar = c(5,5,4,2) + 0.1)
  plot(0, 0, type = "n", xlim = c(0, max(lambda)), ylim = c(0,1),
       xlab = expression(lambda),
       ylab = expression(alpha), cex.lab = 2,
       main = "Entire Regularization Path for ROC-SVM")
  for (i in 1:N) lines(lambda, alpha[i,], type = "l", ...)
  abline(v = lambda, col = "gray", lty=2)
}
