source("r_Gaussian_mixed.R")
N <- 10000 #number of samples
K <- 3 #number of mixed components
p <- 6 #dimension of covariate
mu <- matrix(c(2,2,0,0,0,0,0,0,2,2,0,0,
  0,0,0,0,2,2),nrow = 6)
sigma <- array(c(diag(x=1,nrow = 6),diag(x=1,nrow = 6),diag(x=1,nrow = 6)),dim = c(6,6,3))
weight <- rep(1,3)/3 
data <- r_Gaussian_mixed(n = N, mu = mu, sigma = sigma, weight = weight)


##EM algorithm
source("EM_alg.R")
#EM algorithm without the assumption that the covariance matrix is identity
epsilon <- 1e-3 
t <- Sys.time()
EM_alg(data = data, K = K, epsilon = epsilon)
Sys.time()-t

#EM algorithm with the assumption that the covariance matrix is identity
epsilon <- 1e-3
t <- Sys.time()
EM_alg_identity(data = data, K = K, epsilon = epsilon)
Sys.time()-t

##LearnGMM
source("LearnGMM.R")
delta <- 1e-3
t <- Sys.time()
LearnGMM(data = data, K = K, delta = delta)
Sys.time()-t