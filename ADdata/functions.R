###Broken stick with linear SD transformation
bMSD <- function(x){
  (nrow <- dim(x)[1])
  (ncol <- dim(x)[2])
  (y <- matrix(nrow=nrow,ncol=ncol, dimnames=list(rownames(x), colnames(x))))
  (total <- 100)
  (y[1,2] <- (x[1,2])/100)
  (y[1,1] <- (x[1,1])/100)
  for(i in 2:nrow){
    (y[i,2] <- x[i,2]/(total - x[i-1,1]))
    (total <- total - x[i-1,1])
    (y[i,1] <- x[i,1]/total)
  }
  return(y)
}

###Beta distribution
getBetaParams <- function(mean, sd){
  v <- sd^2
  m <- mean * (1-mean)
  n <- (m/v) - 1
  alpha <- mean*n
  beta <- (1 - mean)*n 
  params <- list(type=1, a=alpha, b=beta, location=0, scale=1)
  return(params)
}