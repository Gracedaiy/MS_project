# Used library
library(mvtnorm)
library(expm)
library(MASS)
library(tidyverse)
set.seed(120733)


## ANN function
rrr <- function(X,Y,lambda,gamma = 2){
  rrresult <- list()
  
  #1 estimate C using LS estimator (generalized inverse is used)
  Clhat <- ginv(t(X)%*%X)%*%t(X)%*% Y
  
  #2 eigendecomposite of CLhar
  Uhat <- svd(X%*%Clhat)$u
  vhat <- svd(X%*%Clhat)$v
  D <- diag(svd(X%*%Clhat)$d)
  
  #3. weight for soft-thredshold and rank
  weight <-  (svd(X%*%Clhat)$d)^(-gamma)
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
n <- 1000

# number of metabolites
m <- 100

# number of taxa
k <- 40

# rank of matrix B
r <- 10

# Generate X from Multivariate Normal Distribution (Original method)
rho <- 0.4
Sigma <- diag(1, n)
for (i in 1:(n-1)){
  for (j in (i+1):n){
    Sigma[i, j] <- Sigma[j, i] <- rho^abs(i-j)
  }
}
Sigma.eigen <- eigen(Sigma)
X.raw <- matrix(rnorm(m*n), nrow = m, ncol = n)
# X <-as.matrix(rbind(rep(1, n), X.raw %*% diag(sqrt(pmax(Sigma.eigen$values, 0))) %*%
#                      t(Sigma.eigen$vectors)))
X <-X.raw %*% diag(sqrt(pmax(Sigma.eigen$values, 0))) %*%
                      t(Sigma.eigen$vectors)

# Generate X from Multivariate Noromal with identity variance and covariance matrix
# Method 1:
# X <- matrix(rnorm(m*n), nrow = m, ncol = n)

# Each entry is iid distributed from N(0, 1), center it afterwards
# Method 2: (don't know why it does not work!)
# X.raw <- matrix(rnorm(m*n), nrow = m, ncol = n)
# X.mean <- matrix(rep(rowMeans(X.raw), each = n), ncol = n, byrow = T)
# X <- X.raw - X.mean


# Method 3: covariance structure & centered
# rho <- 0.4
# Sigma <- diag(1, n)
# for (i in 1:(n-1)){
#   for (j in (i+1):n){
#     Sigma[i, j] <- Sigma[j, i] <- rho^abs(i-j)
#   }
# }
# Sigma.eigen <- eigen(Sigma)
# X.raw <- matrix(rnorm(m*n), nrow = m, ncol = n)
# # X <-as.matrix(rbind(rep(1, n), X.raw %*% diag(sqrt(pmax(Sigma.eigen$values, 0))) %*%
# #                      t(Sigma.eigen$vectors)))
# X.raw <-X.raw %*% diag(sqrt(pmax(Sigma.eigen$values, 0))) %*%
#                       t(Sigma.eigen$vectors)
# X.mean <-  matrix(rep(rowMeans(X.raw), each = n), ncol = n, byrow = T)
# X <- X.raw - X.mean


# Method 4: randomize the variance and covaraince matrix
# Sigma <- diag(runif(n), n)
# for (i in 1:(n-1)){
#   for (j in (i+1):n){
#     Sigma[i, j] <- Sigma[j, i] <- runif(1)
#   }
# }
# Sigma.eigen <- eigen(Sigma)
# X.raw <- matrix(rnorm(m*n), nrow = m, ncol = n)
# X <- X.raw %*% diag(sqrt(pmax(Sigma.eigen$values, 0))) %*%
#                       t(Sigma.eigen$vectors)

# Method 5: randomize the variance and covariance matrix, then centered
# Sigma <- diag(runif(n), n)
# for (i in 1:(n-1)){
#   for (j in (i+1):n){
#     Sigma[i, j] <- Sigma[j, i] <- runif(1)
#   }
# }
# Sigma.eigen <- eigen(Sigma)
# X.raw <- matrix(rnorm(m*n), nrow = m, ncol = n)
# X.raw <- X.raw %*% diag(sqrt(pmax(Sigma.eigen$values, 0))) %*%
#                       t(Sigma.eigen$vectors)
# X.mean <-  matrix(rep(rowMeans(X.raw), each = n), ncol = n, byrow = T)
# X <- X.raw - X.mean


# Method 6: Centered use the grand mean
# rho <- 0.4
# Sigma <- diag(1, n)
# for (i in 1:(n-1)){
#   for (j in (i+1):n){
#     Sigma[i, j] <- Sigma[j, i] <- rho^abs(i-j)
#   }
# }
# Sigma.eigen <- eigen(Sigma)
# X.raw <- matrix(rnorm(m*n), nrow = m, ncol = n)
# X.raw <-X.raw %*% diag(sqrt(pmax(Sigma.eigen$values, 0))) %*%
#                       t(Sigma.eigen$vectors)
# grand.mean <- mean(X.raw)
# X.mean <-  matrix(grand.mean, ncol = n, nrow = m)
# X <- X.raw - X.mean


# Generate low-rank matrix B
sing.vals <- 5*((r:1)/r + 1/r)
B <- rortho(k)[,1:r] %*% diag(sing.vals) %*% rortho(m)[1:r,]

# Generate T vector from Unif(0, 10)
T <- matrix(runif(n, min = 0, max = 10), nrow = 1)

# Generate random error
E <- matrix(rnorm(k*n), nrow = k, ncol = n)

# Calculate our new Z = XB + 1T + E
Z <- B %*% X + matrix(rep(1, k), ncol = 1) %*% T + E

# convergence criterion
vareps <- 10^(-3)

# range of lambda
lambda.seq <- exp(seq(9, 10, 1))
rank.est <- rep(0, length(lambda.seq))
converge.iter <- rep(0, length(lambda.seq))

for (i in 1: length(lambda.seq)) {
  # # initial estimated B
  # qrX <- qr(t(X), tol = 10^-7)
  # B_iter <- qr.coef(qrX, t(Z))
  # 
  # # B_iter <- t(B)
  # 
  # # initial estimated T
  # T_iter <- matrix(rowMeans(t(Z) - t(X) %*% B_iter), nrow = 1)
  
  disT <- 1
  disB <- 1
  disXB <- 1
  iter.count <- 1
  
  T_iter <- matrix(rnorm(n, mean = 20,  sd = 8), nrow = 1)
  T_iter <- T_iter - mean(T_iter)
  # T_iter <- T - mean(T)
  Z_iter <- t(Z) - t(T_iter) %*% t(matrix(rep(1, k), ncol = 1))
  B_iter <- rrr(t(X), Z_iter, lambda = lambda.seq[i])$CShat
    
  while ((((disB > vareps) || (disT > vareps))  && ((disXB > vareps) || (disT > vareps))) && (iter.count < 3000)) {
    T_iter_next <- matrix(rowMeans(t(Z) - t(X) %*% B_iter), nrow = 1)
    # T_iter_next <- T
    
    Z_iter <- t(Z) - t(T_iter_next) %*% t(matrix(rep(1, k), ncol = 1))
    B_iter_next <- rrr(t(X), Z_iter, lambda = lambda.seq[i])$CShat
    # B_iter_next <- t(B)
    B_iter_next_rank <- rrr(t(X), Z_iter, lambda = lambda.seq[i])$r

    disT <- sum((T_iter - T_iter_next)^2)
    disB <- sum((B_iter - B_iter_next)^2)
    disXB <- sum((t(B_iter) %*% X - t(B_iter_next) %*% X)^2)
    
    Z.est <-  t(B_iter) %*% X + matrix(rep(1, k), ncol = 1) %*% T_iter
    B_iter <- B_iter_next
    T_iter <- T_iter_next
    iter.count <- iter.count + 1
    # print(iter.count)
    # print(T_iter[1:10])
  }
  rank.est[i] <- B_iter_next_rank
  converge.iter[i] <- iter.count
  
  # save image for B matrix
  png(paste0("plot-B-RCT", i, ".png"))
  plot(B, B_iter)
  dev.off()

  # save image for T matrix
  png(paste0("plot-T-RCT", i, ".png"))
  plot(T, T_iter)
  abline(a = 0, b = 1)
  dev.off()

  # save image for XB matrix
  png(paste0("plot-XB-RCT", i, ".png"))
  plot(B %*% X, t(B_iter) %*% X)
  abline(a = 0, b = 1)
  dev.off()

  # save image for Z matrix
  png(paste0("plot-Z-RCT", i, ".png"))
  plot(Z, Z.est)
  abline(a = 0, b = 1)
  dev.off()
}


#plot(B%*%X, t(B_iter) %*% X)

