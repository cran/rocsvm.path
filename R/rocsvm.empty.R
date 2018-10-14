rocsvm.empty <- function(K, Left, rho)
{
  Kleft <- K[Left,Left,drop = F]
  nl <- length(Left)
  beta0 <- Kleft %*% rep(1, nl)
  index0 <- which.max(beta0)
  lambda0 <- beta0[index0]/rho

  Elbow <- Left[index0]
  Left <- Left[-index0]

  obj <- list(lambda = lambda0,
              Elbow = Elbow,
              Left = Left)
}
