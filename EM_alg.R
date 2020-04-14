library(MASS)
dmvnorm <- function(x, mu, sigma){
  d <- length(x)
  exp(-t(x-mu)%*%solve(sigma)%*%(x-mu)/2)/((2*pi)^(d/2)*sqrt(abs(det(sigma))))
}

EM_alg <- function(data,K,epsilon){
  #initialization
  N <- nrow(data)
  d <- ncol(data)
  index <- sample(N,size = K,replace = FALSE)
  mu <- t(data[index,])
  sigma <- array(rep(cov(data),K),dim = c(d,d,K))
  alpha <- rep(1,K)/K
  p <- matrix(nrow = N, ncol = K)
  for (i in 1:N) {
    for (j in 1:K) {
      p[i,j] <- dmvnorm(x = data[i, ], mu = mu[, j], sigma = sigma[, , j])
    }
  }
  log_l <- sum(log(p%*%alpha))
  
  repeat{
    #E-Step
    w <- t(apply(data.frame(1:N), 1, function(x){
      p[x,]*alpha/c(p[x,]%*%alpha)
    }))
    
    #M-Step
    N_k <- colSums(w)
    alpha <- N_k/N
    mu <- apply(matrix(1:K), 1, function(x){
      w[, x]%*%data/N_k[x]
    })
    weighted_matrix <- function(x){
      sum <- matrix(0,nrow = d, ncol = d)
      for (i in 1:N) {
        sum <- sum + w[i, x]*(data[i, ]-mu[, x])%*%t(data[i, ]-mu[, x])
      }
      return(sum/N_k[x])
    }
    sigma <- array(apply(matrix(1:K),1,
      weighted_matrix),dim = c(d,d,K))
    
    for (i in 1:N) {
      for (j in 1:K) {
        p[i,j] <- dmvnorm(x = data[i, ], mu = mu[, j], sigma = sigma[, , j])
      }
    }
    log_l_new <- sum(log(p%*%alpha))
    if(abs(log_l_new-log_l)<epsilon){break}
    log_l <- log_l_new
  }
  return(list(mu=mu, sigma=sigma, weight=alpha))
}



dmvnorm_identity <- function(x, mu, sigma_sq){
  d <- length(x)
  exp(-t(x-mu)%*%(x-mu)/(2*sigma_sq))/((2*pi*sigma_sq)^(d/2))
}

EM_alg_identity <- function(data,K,epsilon){
  #initialization
  N <- nrow(data)
  d <- ncol(data)
  index <- sample(N,size = K,replace = FALSE)
  mu <- t(data[index,])
  sigma_sq <- rep(mean(apply(data, 1, function(x){t(x-colMeans(data))%*%(x-colMeans(data))/d})),K)
  alpha <- rep(1,K)/K
  p <- matrix(nrow = N, ncol = K)
  for (i in 1:N) {
    for (j in 1:K) {
      p[i,j] <- dmvnorm_identity(x = data[i, ], mu = mu[, j], sigma_sq = sigma_sq[j])
    }
  }
  log_l <- sum(log(p%*%alpha))
  
  repeat{
    #E-Step
    w <- t(apply(data.frame(1:N), 1, function(x){
      p[x,]*alpha/c(p[x,]%*%alpha)
    }))
    
    #M-Step
    N_k <- colSums(w)
    alpha <- N_k/N
    mu <- apply(matrix(1:K), 1, function(x){
      w[, x]%*%data/N_k[x]
    })
    weighted_matrix_identity <- function(x){
      sum <- 0
      for (i in 1:N) {
        sum <- sum + w[i, x]*t(data[i, ]-mu[, x])%*%(data[i, ]-mu[, x])/d
      }
      return(sum/N_k[x])
    }
    sigma_sq <- apply(matrix(1:K),1,weighted_matrix_identity)
    
    for (i in 1:N) {
      for (j in 1:K) {
        p[i,j] <- dmvnorm_identity(x = data[i, ], mu = mu[, j], sigma_sq = sigma_sq[j])
      }
    }
    log_l_new <- sum(log(p%*%alpha))
    if(abs(log_l_new-log_l)<epsilon){break}
    log_l <- log_l_new
  }
  return(list(mu=mu, sigma_sq=sigma_sq, weight=alpha))
}