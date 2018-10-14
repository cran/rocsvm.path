#' Finding an intercept fixed sensitivity or specificity for ROC-SVM
#'
#' Computes an intercept at a specific sensitivity or specificity level from the ROC-SVM model.
#'
#' @param obj The rocsvm.path object
#' @param lambda The regularization parameter that users want in ROC-SVM model.
#' @param sensitivity Sensitivity in ROC curve, which means True Positive Rate (TPR).
#' @param specificity Specificity in ROC curve, which means True Negative Rate (TNR) = 1-FPR.
#' @author Seung Jun Shin, Do Hyun Kim
#' @seealso \code{\link{rocsvm.path}}
#' @examples
#' # library(rocsvm.path)
#' # The 'obj' comes from an example description of rocsvm.path()
#' \donttest{
#' rocsvm.intercept(obj, lambda = 1, sensitivity = 0.9, specificity = 0.1) }
#' @export
#' @importFrom svmpath poly.kernel radial.kernel
rocsvm.intercept <- function(obj, lambda = 1, sensitivity = 0.5, specificity = 0.5){

  x <- obj$x
  y <- obj$y
  kernel <- obj$kernel
  param.kernel <- obj$param.kernel
  prop <- obj$prop

  alpha0 <- rocsvm.get.solution(obj, lambda)
  K  <- rocsvm.compute.K(x, y, x, kernel, param.kernel, prop)
  fitted <-  (K %*% alpha0) / lambda

  out <- rocsvm.curve(fitted, y, new = FALSE)
  tpr.out <- out$tpr ; fpr.out <- out$fpr

  sort.x <- sort(fitted)

  if(sum(sensitivity == tpr.out) != 0){
    intercept.sens <- sort.x[which(sensitivity == tpr.out)]
    intercept.sens <- mean(intercept.sens)
  }else{
    i1 <- max(which(tpr.out > sensitivity))
    i2 <- min(which(tpr.out < sensitivity))
    y1 <- tpr.out[i1] ; y2 <- tpr.out[i2]
    x1 <- sort.x[i1]; x2 <- sort.x[i2]

    slope <- (y2 - y1)/(x2 - x1)
    intercept.sens <- (sensitivity - y1) / slope + x1
  }

  fpr <- 1 - specificity
  if(sum(fpr == fpr.out) != 0){
    intercept.spec <- sort.x[which(fpr == fpr.out)]
    intercept.spec <- mean(intercept.spec)
  }else{
    i1 <- max(which(fpr.out > fpr))
    i2 <- min(which(fpr.out < fpr))
    y1 <- fpr.out[i1] ; y2 <- fpr.out[i2]
    x1 <- sort.x[i1]; x2 <- sort.x[i2]

    slope <- (y2 - y1)/(x2 - x1)
    intercept.spec <- (fpr - y1) / slope + x1
  }
  return(data.frame(lambda = lambda, sensitivity = sensitivity,'intercept.sensitivity' = intercept.sens,
                    specificity = specificity, 'intercept.specificity' = intercept.spec))
}
