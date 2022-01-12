# Used library
library(mvtnorm)
library(expm)
library(rrpack)
library(matlib)
set.seed(120733)

# sample size n (n person's data have been taken)
n <- 50

# number of metabolites
m <- 100

# number of taxa
k <- 40

# rank of matrix B
r <- 10

# H matrix
one.vec <- matrix(rep(1, k), nrow = k)
H.full <- diag(k) - one.vec %*% solve(t(one.vec) %*% one.vec) %*% t(one.vec)
H <- H.full[-1, ]

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
X <-as.matrix(rbind(rep(1, m), X.raw %*% diag(sqrt(pmax(Sigma.eigen$values, 0))) %*% 
                      t(Sigma.eigen$vectors)))

# function to generate a random orthonormal matrix
rortho <- function(x){
  return(qr.Q(qr(array(runif(x), dim = c(x,x)))))
}

# Generate low-rank matrix B
sing.vals <- 5*((r:1)/r + 1/r)
B <- rortho(k)[,1:r] %*% diag(sing.vals) %*% rortho(m+1)[1:r,]


# Generate random error
E <- matrix(rnorm(k*n), nrow = k, ncol = n)

# Calculate our new Y = HZ = HBX + HE
Y <- H %*% B %*% X + H %*% E
C <- H %*% B

# LS esitmates of variance-covaraince matrix
sig.eps.tilda <- (1/n) * (Y %*% t(Y) - Y %*% t(X) %*% Ginv(X %*% t(X)) %*% X %*% t(Y))
C.tilda <- Y %*% t(X) %*% Ginv(X %*% t(X))

var.yx <- (1/n) * (Y %*% t(X))
var.xx <- (1/n) * (X %*% t(X))
sqrt.eps.tilda <- sqrtm(sig.eps.tilda)

V.hat <- svd(solve(sqrt.eps.tilda) %*% var.yx
             %*% Ginv(var.xx) %*% t(var.yx) 
             %*% solve(sqrt.eps.tilda))$u[,1:r]
C.hat <- sqrt.eps.tilda %*% V.hat %*% t(V.hat) %*% solve(sqrt.eps.tilda) %*% C.tilda
Y.hat <-  C.hat %*% X
