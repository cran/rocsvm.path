#' Compute the kernel matrix for ROC-SVM path 
#' @description
#' Compute the kernel matrix for ROC-SVM path. This function comes from \code{svmpath} package by Trevor Hastie. If you want to know details of this function, refer the \code{svmpath} package.
#' @param x An n x p matrix of features 
#' @param y An m x p matrix of features 
#' @param param.kernel The parameter(s) for the kernel. For this radial kernel, the parameter is known in the fields as "gamma". For the polynomial kernel, it is the "degree"
#' @param ... unused 
#' @export
radial.kernel <- function (x, y = x, param.kernel = 1/p, ...) 
{
  n <- nrow(x)
  m <- nrow(y)
  p <- ncol(x)
  normx <- drop((x^2) %*% rep(1, p))
  normy <- drop((y^2) %*% rep(1, p))
  a <- x %*% t(y)
  a <- (-2 * a + normx) + outer(rep(1, n), normy, "*")
  exp(-a * param.kernel)
}