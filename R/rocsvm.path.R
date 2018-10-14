#' Fit the entire regularization path for ROC-Support Vector Machine (ROC-SVM)
#' @description
#' This algorithm computes the entire regularization path for the ROC-Support Vector Machine with a relatively low cost compared to quadratic programming problem.
#'
#' @param x The data matrix (n x p) with n rows (observations) on p variables (columns)
#' @param y The \code{{-1, 1}} valued response variable.
#' @param kernel This is a user-defined function. Provided options are polynomial kernel; \code{poly} (the default, with parameter set to default to a linear kernel) and radial kernel; \code{radial}.
#' @param param.kernel The parameter(s) for the kernel. For this radial kernel, the parameter is known in the fields as "gamma". For the polynomial kernel, it is the "degree"
#' @param lambda.min The smallest value of lambda for termination of the algorithm (the default is \code{lambda.min = 1e-05}).
#' @param prop The proportion of large class corresponding a point of small class by speed-up tricks (the default is \code{prop = 0.5}). If you don't want to use the "speed-up tricks", then set \code{prop} to 1.
#' @param rho A positive constant
#' @param Nmoves The maximum number of iterations the rocsvm.path algorithm
#' @param eps An adjustment computing errors
#' @return A 'rocsvm.path' object is returned, for which there are \code{lambda} values and corresponding values of \code{alpha} for each data point.
#' @author Seung Jun Shin, Do Hyun Kim
#' @seealso  \code{\link{rocsvm.get.solution}}, \code{\link{plot.rocsvm}},  \code{\link{rocsvm.intercept}}
#' @examples
#' library(rocsvm.path)
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
#' obj <- rocsvm.path(x, y, rho, kernel, param.kernel, prop)
#' @export
rocsvm.path <- function(x, y, rho = 1,
                        kernel = poly.kernel, param.kernel = 1,
                        prop = .5,
                        lambda.min = 1.0e-5,
                        eps = 1.0e-5,
                        Nmoves = 500)
{

p <- ncol(x)
n <- length(y)

x.ps <- x[y == 1,,drop = F]
x.mn <- x[y != 1,,drop = F]

n.ps <- sum(y == 1)
n.mn <- sum(y != 1)


X.ps <- X.mn <- NULL
sub.n <- max(1, floor(prop * n.mn))
id.mat <- matrix(0, n.ps, sub.n)
for (i in 1:n.ps) {
  o.id <- order(apply((t(x.mn) - x.ps[i,])^2, 2, sum))
  id.mat[i,] <- sort(o.id[1:sub.n])
  X.mn <- rbind(X.mn, x.mn[id.mat[i,],])
  X.ps <- rbind(X.ps, x.ps[rep(i, sub.n),])
}

N <- nrow(X.mn)

K.pp <- kernel(X.ps, X.ps, param.kernel)
K.pn <- kernel(X.ps, X.mn, param.kernel)
K.np <- kernel(X.mn, X.ps, param.kernel)
K.nn <- kernel(X.mn, X.mn, param.kernel)

K <- K.pp - K.pn - K.np + K.nn

Elbow.list <- as.list(seq(Nmoves))
alpha <- matrix(0, N, Nmoves)
lambda <- NULL

init <- rocsvm.initialization(K, N, rho)

alpha[,1] <- init$alpha
lambda[1] <- lambda0 <- init$lambda
Elbow.list[[1]] <- Elbow <- init$Elbow
Right <- init$Right
Left <- init$Left

k <- 1
fl <- (K %*% alpha[,1])/lambda0

# update
while ((k < Nmoves) && (lambda[k] > lambda.min)) {
  if (length(Elbow) == 0) { # empty elbow
    temp <- rocsvm.empty(K, Left, rho)
    lambda[k + 1] <- temp$lambda
    Elbow <- temp$Elbow
    Left <- temp$Left
    alpha[, k + 1] <- alpha[, k]
    fl <- (K %*% alpha[,k + 1])/lambda[k + 1]
  } else { #not empty
    ne <- length(Elbow)
    one <- rep(1, ne)
    Kstar <- K[Elbow, Elbow, drop = F]
    K.inv <- solve(Kstar)
    b <- rho * K.inv %*% one
    hl <- K[, Elbow, drop = FALSE] %*% b
    dl <- fl - hl


    # Exit from Elbow
    lambda.left <- lambda[k] - (alpha[Elbow, k] - 1)/b
    lambda.left[abs(b) < eps] <- -1
    lambda.right <- lambda[k] - alpha[Elbow, k]/b
    lambda.right[abs(b) < eps] <- -1
    lambda01 <- c(lambda.right, lambda.left)
    lambda.exit <- max(lambda01[lambda01 < (lambda[k] - eps)], na.rm = T)

    # Enter to Elbow
    lambdai <- lambda[k] * (dl)/(rho - hl)
    lambdai[Elbow] <- -Inf
    lambda.entry <- max(c(lambdai[lambdai > 0 & lambdai < (lambda[k] - eps)], -Inf), na.rm= T)


    # Find next ell
    lambda.max <- max(lambda.entry, lambda.exit)

    # Update values
    lambda[k + 1] <- lambda.max
    alpha[, k + 1] <- alpha[, k]
    alpha[Elbow, k + 1] <- alpha[Elbow, k] - (lambda[k] - lambda[k + 1]) * b
    fl <- (lambda[k]/lambda[k + 1]) * dl + hl

    # Update Sets
    if (lambda.entry > lambda.exit) { # Enter to Elbow
      i <- match(lambda.entry, lambdai, 0)[1]
      if (match(i, Left, FALSE)) {
        Left <- setdiff(Left, i)
      } else {
        Right <- setdiff(Right, i)
      }
      Elbow <- c(Elbow, i)
    } else { # Exit from Elbow
      idrop <- Leaveright <- NULL
      i <- Elbow[abs(lambda.right - lambda.exit) < eps]
      if (length(i) > 0) {
        Leaveright <- rep(TRUE, length(i))
        idrop <- i
      }
      i <- Elbow[abs(lambda.left - lambda.exit) < eps]
      if (length(i) > 0) {
        Leaveright <- c(Leaveright, rep(FALSE, length(i)))
        idrop <- c(idrop, i)
      }
      for (j in seq(along = idrop)) {
        if (Leaveright[j]) {
          Right <- c(Right, idrop[j])
        } else {
          Left <- c(Left, idrop[j])
        }
        mi <- match(idrop[j], Elbow)
        Elbow <- Elbow[-mi]
      }
    }
  }
  k <- k + 1
  Elbow.list[[k]] <- Elbow
  if (min(alpha[,k]) < -1.0e-5) break
}
lambda <- c(lambda[1:(k-1)], 0)
alpha <- cbind(alpha[,1:(k-1)], rep(0, N))
Elbow <- Elbow.list[1:(k-1)]

obj <- list(lambda = lambda,
            alpha = alpha,
            Elbow = Elbow,
            K = K,
            kernel = kernel,
            param.kernel = param.kernel,
            neighbor.mn = id.mat,
            x = x, y = y, prop = prop)
class(obj)<-"rocsvm"
obj
}
