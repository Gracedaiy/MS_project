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
X2 <- X[, -1]

# Generate low-rank matrix B (m+1, k)
# sing.vals <- 5*((r:1)/r + 1/r)
# B2 <- t(rortho(k)[,1:r] %*% diag(sing.vals) %*% rortho(m)[1:r,])
C_0 <- matrix(rnorm(m*r), nrow = m)
C_1 <- matrix(rnorm(r*k), nrow = r)
sr <- 0.8
B2 <- sr * (C_0 %*% C_1)

# Generate T vector from Unif(0, 10), (n, 1)
T <- matrix(runif(n, min = 0, max = 10), ncol = 1)
# T <- T - mean(T)

# Generate b vector for intercept from N(0,1), (k, 1)
b <- matrix(rnorm(k), nrow = 1)
B <- as.matrix(rbind(b, B2))

# Generate random error
E <- matrix(rnorm(k*n), nrow = n, ncol = k, byrow = T)

# Calculate our new Z = XB + 1T + E, (n, k)
Z <- X %*% B + matrix(rep(T, each = k), nrow = n, byrow = T) + E

# convergence criterion
vareps <- 10^(-3)

# Validation set1
E_v1 <- matrix(rnorm(k*n), nrow = n, ncol = k, byrow = T)
Z_v1 <- X %*% B + matrix(rep(T, each = k), nrow = n, byrow = T) + E_v1

# Validation set2
E_v2 <- matrix(rnorm(k*n), nrow = n, ncol = k, byrow = T)
Z_v2 <- X %*% B + matrix(rep(T, each = k), nrow = n, byrow = T) + E_v2

# Validation set3
E_v3 <- matrix(rnorm(k*n), nrow = n, ncol = k, byrow = T)
Z_v3 <- X %*% B + matrix(rep(T, each = k), nrow = n, byrow = T) + E_v3

# Validation set4
E_v4 <- matrix(rnorm(k*n), nrow = n, ncol = k, byrow = T)
Z_v4 <- X %*% B + matrix(rep(T, each = k), nrow = n, byrow = T) + E_v4

# Validation set5
E_v5 <- matrix(rnorm(k*n), nrow = n, ncol = k, byrow = T)
Z_v5 <- X %*% B + matrix(rep(T, each = k), nrow = n, byrow = T) + E_v5


# range of lambda
lambda.seq <- exp(seq(1, 4, 0.2))
rank.est <- rep(0, length(lambda.seq))
converge.iter <- rep(0, length(lambda.seq))
num_obj <- rep(0, length(lambda.seq))
pred_err1 <- rep(0, length(lambda.seq))
pred_err2 <- rep(0, length(lambda.seq))
pred_err3 <- rep(0, length(lambda.seq))
pred_err4 <- rep(0, length(lambda.seq))
pred_err5 <- rep(0, length(lambda.seq))


start_time <- Sys.time()

for (i in 1: length(lambda.seq)) {
  # initial estimated B
  # qrX <- qr(X, tol = 10^-7)
  # B_iter <- qr.coef(qrX, Z)
  B_iter <- ginv(t(X)%*%X)%*%t(X)%*% Z
  
  # initial estimated T (RLS)
  T_iter <- matrix(rowMeans(Z - X %*% B_iter), ncol = 1)
  T_iter <- T_iter - mean(T_iter)
  
  disT <- 1
  disB <- 1
  disXB <- 1
  iter.count <- 1
  
  while ((((disB > vareps) || (disT > vareps))  && ((disXB > vareps) || (disT > vareps))) && (iter.count < 3000)) {
    Z_iter_B <- Z - matrix(rep(T_iter, each = k), nrow = n, byrow = T)
    
    B_iter_next <- rrr(X, Z_iter_B, lambda = lambda.seq[i])$CShat
    B_iter_next_rank <- rrr(X, Z_iter_B, lambda = lambda.seq[i])$r
    
    # Using RLS method
    T_iter_next <- matrix(rowMeans(Z - X %*% B_iter), ncol = 1)
    T_iter_next <- T_iter_next - mean(T_iter_next)
    
    iter_weight <- sort(rrr(X, Z_iter_B, lambda = lambda.seq[i])$weight)
    disT <- sum((T_iter - T_iter_next)^2)
    disB <- sum((B_iter - B_iter_next)^2)
    disXB <- sum((X %*% B_iter - X %*% B_iter_next)^2)
    # print(c(disT, disB, disXB))
    # print(sum(T_iter))
    
    B_iter <- B_iter_next
    T_iter <- T_iter_next
    iter.count <- iter.count + 1
    # print(iter.count)
  }
  rank.est[i] <- B_iter_next_rank
  converge.iter[i] <- iter.count
  num_obj[i] <- 1/2 * sum((Z - X %*% B_iter - matrix(rep(T_iter, each = k), nrow = n, byrow = T))^2) + 
    lambda.seq[i] * sum(iter_weight*svd(X %*% B_iter)$d)
  Z.est <-  X %*% B_iter + matrix(rep(T_iter, each = k), nrow = n, byrow = T)
  pred_err1[i] <- sum((Z_v1 - Z.est)^2)
  pred_err2[i] <- sum((Z_v2 - Z.est)^2)
  pred_err3[i] <- sum((Z_v3 - Z.est)^2)
  pred_err4[i] <- sum((Z_v4 - Z.est)^2)
  pred_err5[i] <- sum((Z_v5 - Z.est)^2)
 
  # # save image for B matrix
  # png(paste0("plot-B-rls2-wocT-newB.png"))
  # plot(B2, B_iter[-1, ])
  # abline(a = 0, b = 1)
  # dev.off()
  # 
  # # save image for XB matrix
  # png(paste0("plot-XB-rls2-wocT-newB.png"))
  # plot(X2 %*% B2,  X2 %*% B_iter[-1,])
  # abline(a = 0, b = 1)
  # dev.off()
  # 
  # # save image for T matrix
  # png(paste0("plot-T-rls2-wocT-newB.png"))
  # plot(T, T_iter)
  # abline(a = 0, b = 1)
  # dev.off()
  # 
  # # save image for T matrix
  # png(paste0("plot-beta-rls2-wocT-newB.png"))
  # plot(b, B_iter[1,])
  # abline(a = 0, b = 1)
  # dev.off()
  # 
  # # save image for Z matrix
  # png(paste0("plot-Z-rls2-wocT-newB.png"))
  # plot(Z, Z.est)
  # abline(a = 0, b = 1)
  # dev.off()
}
pred_err <- 1/5 * (pred_err1 + pred_err2 + pred_err3 + pred_err4 + pred_err5)
end_time <- Sys.time()
run_time <- end_time - start_time

T[1:10]
T_iter[1:10]
b[1:10]
B_iter[1,][1:10]
B2[1:2, 1:7]
(B_iter[-1,])[1:2, 1:7]
(X2 %*% B2)[1:2, 1:7]
(X2 %*% B_iter[-1,])[1:2, 1:7]