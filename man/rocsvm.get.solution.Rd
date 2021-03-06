% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rocsvm.get.solution.R
\name{rocsvm.get.solution}
\alias{rocsvm.get.solution}
\title{Finding solutions fixed the regularization parameter of ROC-SVM.}
\usage{
rocsvm.get.solution(obj, lambda)
}
\arguments{
\item{obj}{The rocsvm.path object}

\item{lambda}{The regularization parameter that users want in ROC-SVM model.}
}
\description{
Computes solution \code{alpha} values from a fixed regularization parameter, \code{lambda} value for ROC-SVM path object.
}
\examples{
# library(rocsvm.path)
# The 'obj' comes from an example description of rocsvm.path()
\donttest{
rocsvm.get.solution(obj, lambda = 1)
}

}
\seealso{
\code{\link{rocsvm.path}}
}
\author{
Seung Jun Shin, Do Hyun Kim
}
