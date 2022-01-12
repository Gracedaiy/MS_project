# Used library
library(mvtnorm)
library(expm)
library(MASS)
library(tidyverse)
set.seed(120733)
# set.seed(19911130)

## ANN function
rrr <- function(X,Y,lambda,gamma = 0){
  rrresult <- list()
  
  #1 estimate C using LS estimator (generalized inverse is used)
  Clhat <- ginv(t(X)%*%X)%*%t(X)%*% Y
  
  #2 eigendecomposite of CLhar
  Uhat <- svd(X%*%Clhat)$u
  vhat <- svd(X%*%Clhat)$v
  D <- diag(svd(X%*%Clhat)$d)
  
  #3. weight for soft-thredshold and rank
  weight <-  (svd(X%*%Clhat)$d)^(-gamma)
  rrresult$weight <- weight
  rrresult$r <-  sum(svd(X%*%Clhat)$d > lambda^(1/(gamma+1)))
  
  #4. soft-thredsholding
  diagvalue <-  svd(X%*%Clhat)$d -lambda*weight
  th <- ifelse(diagvalue < 0, 0, diagvalue)
  SDhat <-  diag(th)
  
  rrresult$thd <- diagvalue
  
  #5. low rank solution C
  rrresult$CShat <-  Clhat %*% vhat %*% ginv(D) %*% SDhat %*% t(vhat)
  
  return(rrresult)
}

# function to generate a random orthonormal matrix
rortho <- function(x){
  return(qr.Q(qr(array(runif(x), dim = c(x,x)))))
}

# simulation settings
# sample size n (n person's data have been taken)
n <- 100

# number of metabolites
m <- 100

# number of taxa
k <- 40

# rank of matrix B
r <- 10

# Generate X from Multivariate Noromal with identity variance and covariance matrix (n, m+1)
X <- matrix(c(rep(1, n), rnorm(m*n)), nrow = n, ncol = m+1)
X2 <- X[,-1]

# Generate low-rank matrix B (m+1, k)
# sing.vals <- 5*((r:1)/r + 1/r)
# B2 <- t(rortho(k)[,1:r] %*% diag(sing.vals) %*% rortho(m)[1:r,])
C_0 <- matrix(rnorm(m*r), nrow = m)
C_1 <- matrix(rnorm(r*k), nrow = r)
sr <- 0.4
B2 <- sr * (C_0 %*% C_1)

# Generate T vector from Unif(0, 10), (n, 1)
T <- matrix(runif(n, min = 0, max = 10), ncol = 1)
# T <- T - mean(T)

# Generate b vector for intercept from N(0,1), (k, 1)
b <- matrix(rnorm(k), nrow = 1)
B <- as.matrix(rbind(b, B2))

# Generate random error
E <- matrix(rnorm(k*n), nrow = n, ncol = k, byrow = T)

# Calculate our new Z = XB + E, (n, k)
Z <- X %*% B + matrix(rep(T, each = k), nrow = n, byrow = T)+ E

# convergence criterion
vareps <- 10^(-3)

# Validation set
E_v <- matrix(rnorm(k*n), nrow = n, ncol = k, byrow = T)
Z_v <- X %*% B + matrix(rep(T, each = k), nrow = n, byrow = T) + E_v

# range of lambda
lambda.seq <- exp(seq(1, 4, 0.2))
rank.est <- rep(0, length(lambda.seq))
num_obj <- rep(0, length(lambda.seq))
pred_err <- rep(0, length(lambda.seq))

start_time <- Sys.time()

for (i in 1: length(lambda.seq)) {
  B_est <- rrr(X, Z, lambda = lambda.seq[i])$CShat
  B_est_rank <- rrr(X, Z, lambda = lambda.seq[i])$r
  iter_weight <- sort(rrr(X, Z, lambda = lambda.seq[i])$weight)
  
  rank.est[i] <- B_est_rank
  num_obj[i] <- 1/2 * sum((Z - X %*% B_est)^2) + lambda.seq[i] * sum(iter_weight*svd(X %*% B_est)$d)
  Z.est <-  X %*% B_est
  pred_err[i] <- sum((Z_v - Z.est)^2)
  
  # # save image for B matrix
  # png(paste0("plot-B-rrr-wocT-newB.png"))
  # plot(B2, B_est[-1, ])
  # abline(a = 0, b = 1)
  # dev.off()
  # 
  # # save image for XB matrix
  # png(paste0("plot-XB-rrr-wocT-newB.png"))
  # plot(X2 %*% B2,  X2 %*% B_est[-1,])
  # abline(a = 0, b = 1)
  # dev.off()
  # 
  # # # save image for T matrix
  # # png(paste0("plot-T-ls2-cT.png"))
  # # plot(T, T_iter)
  # # abline(a = 0, b = 1)
  # # dev.off()
  # 
  # # save image for T matrix
  # png(paste0("plot-beta-rrr-wocT-newB.png"))
  # plot(b, B_est[1,])
  # abline(a = 0, b = 1)
  # dev.off()
  # 
  # # save image for Z matrix
  # png(paste0("plot-Z-rrr-wocT-newB.png"))
  # plot(Z, Z.est)
  # abline(a = 0, b = 1)
  # dev.off()
}

end_time <- Sys.time()
run_time <- end_time - start_time


b[1:10]
B_est[1,][1:10]
B2[1:2, 1:7]
(B_est[-1,])[1:2, 1:7]
(X2 %*% B2)[1:2, 1:7]
(X2 %*% B_est[-1,])[1:2, 1:7]