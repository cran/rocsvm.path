rocsvm.curve <- function(x , y, new = T){

  sort.x <- sort(x)
  n <- length(y)
  pos.n <- length(which(y == 1))
  neg.n <- length(which(y != 1))
  pred <- matrix(0, ncol = n, nrow = n)
  tpr <- fpr <- rep(0, n)

  for(i in 1:n){
    pred[,i] <- ifelse(x >= sort.x[i], 1, min(y))
    if(i == n){
      pred[, i] <- rep(min(y), n)
    }
  }

  for(i in 1:n){
    tpr[i] <- length(which(y == pred[,i] & y == 1)) / pos.n
    fpr[i] <- length(which(y != pred[,i] & y != 1)) / neg.n
  }

  if(new == T){
    plot(fpr, tpr, type= 'l', ylim =c(0,1), xlim=c(0,1))
    abline(0, 1, lty =2, col ='lightgray')
  }
  return(list(tpr = tpr, fpr= fpr))
}


