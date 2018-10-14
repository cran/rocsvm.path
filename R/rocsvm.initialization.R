rocsvm.initialization <- function(K, N, rho)
{
  beta0 <- K %*% rep(1, N)
  index0 <- which.max(beta0)
  lambda0 <- beta0[index0]/rho

  Left <- setdiff(1:N, index0)
  Right <- NULL
  Elbow <- index0

  obj <- list(lambda = lambda0,
              alpha = rep(1, N),
              Elbow = Elbow,
              Left = Left,
              Right = Right)
}
