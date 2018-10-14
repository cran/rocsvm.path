rocsvm.compute.K <- function(x, y, new.x = x, kernel = poly.kernel, param.kernel = 1, prop = 1)
{
  p <- ncol(x)
  n <- length(y)

  if (length(new.x) == 1) new.x <- matrix(new.x, p, 1)

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

  K.p <- kernel(new.x, X.ps, param.kernel)
  K.n <- kernel(new.x, X.mn, param.kernel)

  return(K.p - K.n)
}
