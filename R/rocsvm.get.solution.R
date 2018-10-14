#' Finding solutions fixed the regularization parameter of ROC-SVM.
#'
#' Computes solution \code{alpha} values from a fixed regularization parameter, \code{lambda} value for ROC-SVM path object.
#'
#' @param obj The rocsvm.path object
#' @param lambda The regularization parameter that users want in ROC-SVM model.
#' @author Seung Jun Shin, Do Hyun Kim
#' @seealso \code{\link{rocsvm.path}}
#' @examples
#' # library(rocsvm.path)
#' # The 'obj' comes from an example description of rocsvm.path()
#' \donttest{
#' rocsvm.get.solution(obj, lambda = 1)
#' }
#'
#' @export
rocsvm.get.solution <- function(obj, lambda)
{
  lambda.grid <- obj$lambda
  alpha  <- obj$alpha
  N <- nrow(alpha)
  if (lambda >= lambda.grid[1]) {
    solution <- rep(1, N)
  } else if (lambda < 0) {
    stop("lambda should be positive!")
  } else {
    id1 <- sum(lambda.grid > lambda)
    id2 <- id1 + 1

    lambda1 <- lambda.grid[id1]
    lambda2 <- lambda.grid[id2]

    alpha1 <- alpha[,id1]
    alpha2 <- alpha[,id2]

    slope <- (alpha1 - alpha2)/(lambda1 - lambda2)
    solution <- alpha2 + slope * (lambda - lambda2)
  }
  return(solution)
}
