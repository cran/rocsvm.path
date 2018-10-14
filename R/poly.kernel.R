#' Compute the kernel matrix for ROC-SVM path 
#' @description
#' Compute the kernel matrix for ROC-SVM path. This function comes from \code{svmpath} package by Trevor Hastie. If you want to know details of this function, refer the \code{svmpath} package. 
#' @param x An n x p matrix of features 
#' @param y An m x p matrix of features 
#' @param param.kernel The parameter(s) for the kernel. For the radial kernel, the parameter is known in the fields as "gamma". For the polynomial kernel, it is the "degree"
#' @param ... unused 
#' @export
poly.kernel <- function (x, y = x, param.kernel = 1, ...) 
{
  if (is.null(param.kernel)) 
    param.kernel <- 1
  K = if (param.kernel == 1) 
    x %*% t(y)
  else (x %*% t(y) + 1)^param.kernel
  if (param.kernel == 1) 
    attr(K, "linear") = TRUE
  K
}