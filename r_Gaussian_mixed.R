r_Gaussian_mixed <- function(n=100,mu,sigma,weight,tag=FALSE){
  library(MASS)
  r <- runif(n,0,1)
  c.weight <- cumsum(weight)
  if (tag){
    r_G_m <- matrix(0,ncol = nrow(as.matrix(mu))+1,nrow = n)
    i <- 1
    for (x in r){
      order <- sum(c.weight<x)+1
      r_G_m[i,] <- c(mvrnorm(1,mu[,order],sigma[,,order]),order)
      i=i+1
    }
  }else{
    r_G_m <- matrix(0,ncol = nrow(as.matrix(mu)),nrow = n)
    i <- 1
    for (x in r){
      order <- sum(c.weight<x)+1
      r_G_m[i,] <- mvrnorm(1,mu[,order],sigma[,,order])
      i=i+1
    }
  }
  return(r_G_m)
}
