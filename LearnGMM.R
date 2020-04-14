library(MASS)
library(expm)
pick_uniform <- function(d){
  x <- rnorm(d, 0, 1)
  u <- runif(1, 0, 1)
  rho <- u^(1/d)
  return(rho*x/norm(x,"2"))
}

pairbypair_minus <- function(n){
  M <- choose(n,2)
  D <- matrix(0,nrow = M,ncol = n)
  k <- 1
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      D[k,i] <- 1
      D[k,j] <- -1
      k <- k+1
    }
  }
  return(D)
}


LearnGMM <- function(data, K, delta){
  N <- nrow(data)
  d <- ncol(data)
  Index <- sample(N, size = N/2, replace = FALSE)
  data1 <- data[Index, ]
  data2 <- data[-Index, ]
  N1 <- nrow(data1)
  N2 <- nrow(data2)
  mu_hat <- colMeans(data1)
  M2_cal_hat <- matrix(0, nrow = d, ncol = d)
  for (i in 1:N1) {
    M2_cal_hat <-  M2_cal_hat + data1[i, ]%o%data1[i, ]
  }
  M2_cal_hat <-  M2_cal_hat/N1
  sigma_sq_hat <- eigen(M2_cal_hat-mu_hat%o%mu_hat)$values[K]
  svd2use <- svd(M2_cal_hat-sigma_sq_hat*diag(d))
  M2_hat <- matrix(0, nrow = d, ncol = d)
  for (i in 1:K) {
    M2_hat <- M2_hat + svd2use$d[i]*svd2use$u[, i]%o%svd2use$v[, i]
  }
  U_hat <- svd(M2_hat)$u[, 1:K]
  W_hat <- U_hat%*%sqrtm(ginv(t(U_hat)%*%M2_hat%*%U_hat))
  B_hat <- U_hat%*%sqrtm(t(U_hat)%*%M2_hat%*%U_hat)
  
  whitened_mu_hat <- c(t(W_hat)%*%colMeans(data2))
  
  t <- ceiling(log2(1/delta))
  criterion <- 0
  for (j in 1:t) {
    theta <- pick_uniform(K)
    M3_cal_theta_hat <- matrix(0, nrow = K, ncol = K)
    for (i in 1:N2) {
      whitened_x <- c(t(W_hat)%*%data2[i, ])
      M3_cal_theta_hat <- M3_cal_theta_hat + c(theta%*%whitened_x)*whitened_x%o%whitened_x
    }
    M3_cal_theta_hat <- M3_cal_theta_hat/N2
    M3_theta_hat <- M3_cal_theta_hat
    for (i in 1:d) {
      M3_theta_hat <- M3_theta_hat - sigma_sq_hat*(c(theta%*%W_hat[i, ])*whitened_mu_hat%o%W_hat[i, ]+c(theta%*%W_hat[i, ])*W_hat[i, ]%o%whitened_mu_hat+c(theta%*%whitened_mu_hat)*W_hat[i, ]%o%W_hat[i, ])
    }
    eigen2use <- eigen(M3_theta_hat)
    v <- eigen2use$vectors
    lambda <- eigen2use$values
    D <- pairbypair_minus(K)
    if(min(abs(c(D%*%lambda,lambda)))>criterion){
      criterion <- min(abs(c(D%*%lambda,lambda)))
      v_best <- v
      lambda_best <- lambda
      theta_best <- theta
    }
  }
  mu_i_hat <- apply(matrix(1:K), 1, function(x){
    lambda_best[x]/c(theta_best%*%v_best[, x])*B_hat%*%v_best[, x]
  })
  w_hat <- ginv(mu_i_hat)%*%mu_hat
  sigma_i_sq <- apply(matrix(1:K), 1, function(x){
    (ginv(mu_i_hat)%*%colMeans(data2))[x]/w_hat[x]
  })
  w_hat <- w_hat/sum(w_hat)

  return(list(mu=mu_i_hat,w=w_hat,sigma=sigma_sq_hat))
}