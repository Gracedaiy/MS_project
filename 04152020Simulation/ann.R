# Used library
library(mvtnorm)
library(expm)
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

# Generate X from Multivariate Normal Distribution
rho <- 0.4
Sigma <- diag(1, n)
for (i in 1:(n-1)){
  for (j in (i+1):n){
    Sigma[i, j] <- Sigma[j, i] <- rho^abs(i-j)
  }
}
Sigma.eigen <- eigen(Sigma)
X.raw <- matrix(rnorm(m*n), nrow = m, ncol = n)
X <-as.matrix(rbind(rep(1, n), X.raw %*% diag(sqrt(pmax(Sigma.eigen$values, 0))) %*% 
            t(Sigma.eigen$vectors)))


# Generate low-rank matrix B
sing.vals <- 5*((r:1)/r + 1/r)
B <- rortho(k)[,1:r] %*% diag(sing.vals) %*% rortho(m+1)[1:r,]

# Generate T vector from Unif(0, 10)
T <- matrix(runif(n, min = 0, max = 10), nrow = 1)

# Generate random error
E <- matrix(rnorm(k*n), nrow = k, ncol = n)

# Calculate our new Z = XB + 1T + E
Z <- B %*% X + matrix(rep(1, k), ncol = 1) %*% T + E

# initial estimated B
qrX <- qr(t(X), tol = 10^-7)
B_init <- qr.coef(qrX, t(Z))

# initial estimated T
T_iter <- matrix(rowMeans(t(Z) - t(X) %*% B_init), nrow = 1)
# T_iter <- matrix(runif(n, min = 0, max = 10), nrow = 1)

# converge error control
vareps <- 10^(-7)

disT <- 1
disB <- 1

B_iter <- B_init
count <- 1

# 
# B_iter <- rrr(Z.iter, t(X), penaltySVD = "ann",
#               modstr = list(gamma = gamma, nlambda = 1, lambda = lambda.seq[j]))$coef
# disB <- sum((B_init - B_iter)^2)


while (((disB > vareps) || (disT > vareps)) && (count < 100)) {
  Z_iter <- t(Z) - t(T_iter) %*% t(matrix(rep(1, k), ncol = 1))
  B_iter_next <- rrr(Z_iter, t(X), penaltySVD = "ann",
                modstr = list(gamma = gamma, nlambda = 1, lambda = 40))$coef
  #B_iter_next <- t(B)
  # B_iter_rank <- rrr(Z_iter, t(X), penaltySVD = "ann",
  #                    modstr = list(gamma = gamma, nlambda = 1, lambda = 1))$rank
  T_iter_next <- matrix(rowMeans(t(Z) - t(X) %*% B_iter_next), nrow = 1)
  # T_iter_next <- T
  disT <- sum(abs(T_iter - T_iter_next))
  disB <- sum(abs(B_iter - B_iter_next))
  
  B_iter <- B_iter_next
  T_iter <- T_iter_next
  count <- count + 1
  print(count)
  print(T_iter[1:10])
}

Z.est <-  t(B_iter) %*% X + matrix(rep(1, k), ncol = 1) %*% T_iter

#plot(B%*%X, t(B_iter) %*% X)

