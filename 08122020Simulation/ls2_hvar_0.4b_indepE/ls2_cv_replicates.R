# Used library
library(mvtnorm)
library(expm)
library(MASS)
# library(tidyverse)
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

evaluation <- function(k){ 
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
  write.csv(X, file = paste0('X.csv'), row.names = FALSE)
  # Generate low-rank matrix B (m+1, k)
  # sing.vals <- 5*((r:1)/r + 1/r)
  # B2 <- t(rortho(k)[,1:r] %*% diag(sing.vals) %*% rortho(m)[1:r,])
  C_0 <- matrix(rnorm(m*r), nrow = m)
  C_1 <- matrix(rnorm(r*k), nrow = r)
  sr <- 0.4
  B2 <- sr * (C_0 %*% C_1)
  
  # Generate T vector from Unif(0, 40), (n, 1)
  T <- matrix(runif(n, min = 0, max = 40), ncol = 1)
  write.csv(T, file = paste0('TrueT.csv'), row.names = FALSE)
  # T <- T - mean(T)
  
  # Generate b vector for intercept from N(0,1), (k, 1)
  b <- matrix(rnorm(k), nrow = 1)
  B <- as.matrix(rbind(b, B2))
  write.csv(B, file = paste0('TrueB.csv'), row.names = FALSE)
  
  # # Generate random error
  # E <- matrix(rnorm(k*n), nrow = n, ncol = k, byrow = T)
  # 
  # # Calculate our new Z = XB + 1T + E, (n, k)
  # Z <- X %*% B + matrix(rep(T, each = k), nrow = n, byrow = T) + E
  
  # convergence criterion
  vareps <- 10^(-3)
  
  # range of lambda
  repl <- 50
  lambda.seq <- exp(c(2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.6, 2.65, 2.7, 2.8, 2.9, 3.1, 3.3, 3.5, 3.8))
  # lambda.seq <- exp(seq(1, 4, 0.2))
  # rank.est <- rep(0, length(lambda.seq))
  # converge.iter <- rep(0, length(lambda.seq))
  # num_obj <- rep(0, length(lambda.seq))
  # pred_err <- rep(0, length(lambda.seq))
  estT <- matrix(0, nrow = repl, ncol = n)
  pred_err <- matrix(0, nrow = repl, ncol = length(lambda.seq))
  # opt4.sse <- matrix(0, nrow = length(lambda.seq), ncol = 10)
  # opt5.sse <- matrix(0, nrow = length(lambda.seq), ncol = 10)
  
  rank.est <- rep(0, repl)
  converge.iter <- rep(0, repl)
  num_obj <- rep(0, repl)
  err <- rep(0, repl)
  esterr_B <- rep(0, repl)
  esterr_XB <- rep(0, repl)
  
  start_time <- Sys.time()
  
  for (p in 1:repl) {
    E <- matrix(rnorm(k*n), nrow = n, ncol = k, byrow = T)
    Z <- X %*% B + matrix(rep(T, each = k), nrow = n, byrow = T) + E
    # write.csv(E, file = paste0('E', p, '.csv'), row.names = FALSE)
    # write.csv(Z, file = paste0('Z', p, '.csv'), row.names = FALSE)
    for (j in 1:10) {
      # Validation set
      idx <- sample(1:n, n, replace = F)
      idx_val <- idx[(10*j-9) : (10*j)]
      Z_v <- Z[idx_val,]
      X_v <- X[idx_val,]
      
      # Trainning set
      Z_t <- Z[-idx_val, ]
      X_t <- X[-idx_val, ]
      
      for (i in 1: length(lambda.seq)) {
        # initial estimated B
        # qrX <- qr(X, tol = 10^-7)
        # B_iter <- qr.coef(qrX, Z)
        
        B_iter <- ginv(t(X_t)%*%X_t)%*%t(X_t)%*% Z_t
        
        # initial estimated T (LS)
        T_iter <- matrix(rowMeans(Z_t - X_t %*% B_iter), ncol = 1)
        # T_iter <- T_iter1 - mean(T_iter1)
        
        disT <- 1
        disB <- 1
        disXB <- 1
        iter.count <- 1
        
        while ((((disB > vareps) || (disT > vareps))  && ((disXB > vareps) || (disT > vareps))) && (iter.count < 3000)) {
          Z_iter_B <- Z_t - matrix(rep(T_iter, each = k), nrow = (n-10), byrow = T)
          
          B_iter_next <- rrr(X_t, Z_iter_B, lambda = lambda.seq[i])$CShat
          B_iter_next_rank <- rrr(X_t, Z_iter_B, lambda = lambda.seq[i])$r
          
          # Using LS method
          T_iter_next <- matrix(rowMeans(Z_t - X_t %*% B_iter), ncol = 1)
          
          iter_weight <- sort(rrr(X_t, Z_iter_B, lambda = lambda.seq[i])$weight)
          disT <- sum((T_iter - T_iter_next)^2)
          disB <- sum((B_iter - B_iter_next)^2)
          disXB <- sum((X_t %*% B_iter - X_t %*% B_iter_next)^2)
          # print(c(disT, disB, disXB))
          # print(sum(T_iter))
          
          B_iter <- B_iter_next
          T_iter <- T_iter_next
          iter.count <- iter.count + 1
          # print(iter.count)
        }
        # rank.est[i] <- B_iter_next_rank
        # converge.iter[i] <- iter.count
        # num_obj[i] <- 1/2 * sum((Z - X %*% B_iter - matrix(rep(T_iter, each = k), nrow = n, byrow = T))^2) + 
        #   lambda.seq[i] * sum(iter_weight*svd(X %*% B_iter)$d)
        
        # Method 1: Using ls method to get estimate T by plugging in Z_v
        # T_iter <- matrix(rowMeans(Z_v - X_v %*% B_iter), ncol = 1)
        
        # Method 2: Using mean of the estimated T from iterative algorithm
	# print(c(i, j, mean(T_iter)))
        T_iter <- matrix(rep(mean(T_iter), 10), ncol = 1)
        Z.est <-  X_v %*% B_iter + matrix(rep(T_iter, each = k), nrow = 10, byrow = T)
        Z.est.newval <- X  %*% B_iter + matrix(mean(T_iter), nrow = n, ncol = k )
        # opt4.sse[i, j] <- sum((Z_v - Z.est)^2)
        # opt5.sse[i, j] <- sum((Zv - Z.est.newval)^2)
	# print(c(i, j, mean(T_iter), sum((Z_v - Z.est)^2), sum((Zv - Z.est.newval)^2)))
        pred_err[p,i] <- pred_err[p, i] + sum((Z_v - Z.est)^2)
	# write.csv(B_iter, file = paste0('B', i, '-', j, '.csv'), row.names = FALSE)
      }
    }
    # write.csv(opt4.sse, file = paste0('opt4.sse', p, '.csv'), row.names = FALSE)
    
    tune_lam <- lambda.seq[which.min(pred_err[p,])]
    print(c(p,tune_lam))
    
    B_iter <- ginv(t(X)%*%X)%*%t(X)%*% Z
    
    # initial estimated T (LS)
    T_iter <- matrix(rowMeans(Z- X %*% B_iter), ncol = 1)
    
    disT <- 1
    disB <- 1
    disXB <- 1
    iter.count <- 1
    
    while ((((disB > vareps) || (disT > vareps))  && ((disXB > vareps) || (disT > vareps))) && (iter.count < 3000)) {
      Z_iter_B <- Z - matrix(rep(T_iter, each = k), nrow = n, byrow = T)
      
      B_iter_next <- rrr(X, Z_iter_B, lambda = tune_lam)$CShat
      B_iter_next_rank <- rrr(X, Z_iter_B, lambda = tune_lam)$r
      
      # Using LS method
      T_iter_next <- matrix(rowMeans(Z - X %*% B_iter), ncol = 1)
      # T_iter_next <- T_iter_next - mean(T_iter_next)
      
      iter_weight <- sort(rrr(X, Z_iter_B, lambda = tune_lam)$weight)
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
    rank.est[p] <- B_iter_next_rank
    converge.iter[p] <- iter.count
    num_obj[p] <- 1/2 * sum((Z - X %*% B_iter - matrix(rep(T_iter, each = k), nrow = n, byrow = T))^2) + 
      tune_lam * sum(iter_weight*svd(X %*% B_iter)$d)
    Z.est <-  X %*% B_iter + matrix(rep(mean(T_iter), k*n), nrow = n, byrow = T)
    # fit_err[l] <- sum((Z - Z.est)^2)
    esterr_B[p] <- sum((B2-B_iter[-1, ])^2)
    esterr_XB[p] <- sum((X2 %*% B2 - X2 %*% B_iter[-1, ])^2)
    estT[p,] <- T_iter
    write.csv(B_iter, file = paste0('iterB', p, '.csv'), row.names = FALSE)
    # write.csv(T_iter, file = paste0('iterT', l, '.csv'), row.names = FALSE)
  }
  # write.csv(opt4.sse, file = 'opt4sse.csv', row.names = FALSE)
  # write.csv(opt5.sse, file = 'opt5sse.csv', row.names = FALSE)
  
  write.csv(rank.est, file = 'rank_est.csv', row.names = FALSE)
  write.csv(converge.iter, file = 'cong_iter.csv', row.names = FALSE)
  write.csv(pred_err, file = 'pred_err.csv', row.names = FALSE)
  write.csv(num_obj, file = 'num_obj.csv', row.names = FALSE)
  # write.csv(fit_err, file = 'fit_err.csv', row.names = FALSE)
  write.csv(esterr_B, file = 'fit_B.csv', row.names = FALSE)
  write.csv(esterr_XB, file = 'fit_XB.csv', row.names = FALSE)
  write.csv(estT, file = 'estT.csv', row.names = FALSE)
  end_time <- Sys.time()
  run_time <- end_time - start_time
  print(run_time)
}


args = commandArgs(trailingOnly = TRUE)
evaluation()

