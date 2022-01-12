# Used library
library(mvtnorm)
library(expm)
library(MASS)
library(tidyverse)
set.seed(120733)
# set.seed(19911130)

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
n <- 1000

# number of metabolites
m <- 100

# number of taxa
k <- 40

# rank of matrix B
r <- 10

# Generate X from Multivariate Noromal with identity variance and covariance matrix
X <- matrix(c(rep(1, n), rnorm(m*n)), nrow = m+1, ncol = n, byrow = TRUE)

# Generate low-rank matrix B
sing.vals <- 5*((r:1)/r + 1/r)
B <- rortho(k)[,1:r] %*% diag(sing.vals) %*% rortho(m+1)[1:r,]


# Generate T vector from Unif(0, 10)
T <- matrix(runif(n, min = 0, max = 10), nrow = 1)
T <- T - mean(T)

# Generate random error
E <- matrix(rnorm(k*n), nrow = k, ncol = n)

# Calculate our new Z = XB + 1T + E
Z <- B %*% X + matrix(rep(1, k), ncol = 1) %*% T + E

# convergence criterion
vareps <- 10^(-3)

# range of lambda
lambda.seq <- exp(seq(8, 8, 1))
rank.est <- rep(0, length(lambda.seq))
converge.iter <- rep(0, length(lambda.seq))
num_obj <- rep(0, length(lambda.seq))

# set.seed(Sys.time())
for (i in 1: length(lambda.seq)) {
  # initial estimated B
  
  qrX <- qr(t(X), tol = 10^-7)
  
  # Center the true B
  # B_mean <- matrix(rep(rowMeans(t(B)), each = k), ncol = k, byrow = T)
  # B_iter <- t(B) - B_mean
  
  # randomly generalized starting matrix of B
  B_iter <- qr.coef(qrX, t(Z))
  # B_iter <- B_iter + matrix(rnorm(k * m), ncol = k)
  # B_iter_mean <- matrix(rep(rowMeans(B_iter), each = k), ncol = k, byrow = T)
  # B_iter <- B_iter - B_iter_mean
  
  # initial estimated T
  # T_iter <- matrix(rowMeans(t(Z) - t(X) %*% B_iter), nrow = 1)
  y <- as.vector(Z - t(B_iter) %*% X)
  X_r <- kronecker(diag(n), matrix(rep(1, k), ncol = 1))
  # T_ols <- ginv(t(X_r) %*% X_r) %*% t(X_r) %*% y
  T_ols <- 1/k *(t(X_r) %*% y)
  R <- matrix(rep(1, n), nrow = 1)
  r <- 0
  # T_iter <- T_ols - ginv(t(X_r) %*% X_r) %*% t(R) %*% ginv(R %*% ginv(t(X_r) %*% X_r) %*% t(R)) %*% (R %*% T_ols -r)
  T_iter <- (diag(n) - 1/n*matrix(rep(1, n*n), nrow = n)) %*% T_ols
  
  disT <- 1
  disB <- 1
  disXB <- 1
  iter.count <- 1
  
  while ((((disB > vareps) || (disT > vareps))  && ((disXB > vareps) || (disT > vareps))) && (iter.count < 3000)) {
    Z_iter <- t(Z) - T_iter %*% t(matrix(rep(1, k), ncol = 1))
    B_iter_next <- rrr(t(X), Z_iter, lambda = lambda.seq[i])$CShat
    B_iter_next_rank <- rrr(t(X), Z_iter, lambda = lambda.seq[i])$r
    
    # T_iter_next <- matrix(rowMeans(t(Z) - t(X) %*% B_iter_next), nrow = 1)
    y_iter <- as.vector(Z - t(B_iter) %*% X)
    # T_ols_iter <- ginv(t(X_r) %*% X_r) %*% t(X_r) %*% y_iter
    # T_iter_next <- T_ols_iter - ginv(t(X_r) %*% X_r) %*% t(R) %*% ginv(R %*% ginv(t(X_r) %*% X_r) %*% t(R)) %*% (R %*% T_ols_iter -r)
    T_ols_iter <- 1/k *(t(X_r) %*% y_iter)
    T_iter_next <- (diag(n) - 1/n*matrix(rep(1, n*n), nrow = n)) %*% T_ols_iter
    
    iter_weight <- sort(rrr(t(X), Z_iter, lambda = lambda.seq[i])$weight)
    disT <- sum((T_iter - T_iter_next)^2)
    disB <- sum((B_iter - B_iter_next)^2)
    disXB <- sum((t(B_iter) %*% X - t(B_iter_next) %*% X)^2)

    # Z.est <-  t(B_iter) %*% X + matrix(rep(1, k), ncol = 1) %*% t(T_iter)
    B_iter <- B_iter_next
    T_iter <- T_iter_next
    iter.count <- iter.count + 1
    # print(iter.count)
  }
  rank.est[i] <- B_iter_next_rank
  converge.iter[i] <- iter.count
  num_obj[i] <- 1/2 * sum((t(Z) - t(X) %*% B_iter - T_iter %*% matrix(rep(1, k), nrow = 1))^2) + 
    lambda.seq[i] * sum(iter_weight*svd(t(X) %*% B_iter)$d)
  Z.est <-  t(B_iter) %*% X + matrix(rep(1, k), ncol = 1) %*% t(T_iter)
  # save image for B matrix
  png(paste0("plot-B-rlscT2.png"))
  plot(B, t(B_iter))
  dev.off()

  # save image for T matrix
  png(paste0("plot-T-rlscT2.png"))
  plot(T, T_iter)
  abline(a = 0, b = 1)
  dev.off()

  # save image for XB matrix
  png(paste0("plot-XB-rlscT2.png"))
  plot(B %*% X, t(B_iter) %*% X)
  abline(a = 0, b = 1)
  dev.off()

  # save image for Z matrix
  png(paste0("plot-Z-rlscT2.png"))
  plot(Z, Z.est)
  abline(a = 0, b = 1)
  dev.off()
}
