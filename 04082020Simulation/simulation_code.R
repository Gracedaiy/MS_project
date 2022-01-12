# Used library
library(mvtnorm)
library(expm)
library(rrpack)

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
X <- X.raw %*% diag(sqrt(pmax(Sigma.eigen$values, 0))) %*% 
  t(Sigma.eigen$vectors)

# function to generate a random orthonormal matrix
rortho <- function(x){
  return(qr.Q(qr(array(runif(x), dim = c(x,x)))))
}

# Generate low-rank matrix B
sing.vals <- 5*((r:1)/r + 1/r)
B <- rortho(k)[,1:r] %*% diag(sing.vals) %*% rortho(m)[1:r,]

# Generate T vector from Unif(0, 10)
T <- matrix(runif(n, min = 0, max = 10), nrow = 1)

# Generate random error
E <- matrix(rnorm(k*n), nrow = k, ncol = n)

# Calculate our new Z = XB + 1T + E
Z <- B %*% X + matrix(rep(1, k), ncol = 1) %*% T + E

# Maximum rank you want to use in `rrr` function
maxr <- 20
gamma <- 1

# The user specified upper bound
nrank <- min(k, min(dim(X)), maxr)
rmax <- min(k, min(dim(X)))

# initial estimated B
qrX <- qr(t(X), tol = 10^-7)
B_ls <- qr.coef(qrX, t(Z))

# estimated XC
XB_ls <- qr.fitted(qrX, t(Z))
svdXB_ls <- svd(XB_ls, nu = rmax, nv = rmax)



# svdXB <- svd(t(B %*% X), nu = rmax, nv = rmax)

Ad <- (svdXB_ls$d[1: rmax])^2
Ad <- Ad[Ad > 10^(-7)]
rmax <- length(Ad) 
Adx <- Ad ^ ((1 + gamma) / 2)

# The optimial range for lambda
# lambda.seq <-seq(log(max(Adx)), log(Adx[min(nrank + 1, rmax)]),
#                  length = 100)
# lambda.seq <-exp(seq(log(Adx[2]), log(Adx[min(nrank + 1, rmax)]),
#                                   length = 100))

lambda.seq <- seq(100, 150, 0.5)
rank <- rep(0, length(lambda.seq))
sumsq <- rep(0, length(lambda.seq))

# Run 10 iterations
for (j in 1:length(lambda.seq)) {
  # Starting point
  T.iter <- matrix(rowMeans(t(Z) - t(X) %*% B_ls), nrow = 1)
  for (i in 1:15) {
    Z.iter <- t(Z) - t(T.iter) %*% t(matrix(rep(1, k), ncol = 1))
    B.iter <- rrr(Z.iter, t(X), penaltySVD = "ann",
                  modstr = list(gamma = gamma, nlambda = 1, lambda = lambda.seq[j]))$coef
    B.rank <- rrr(Z.iter, t(X), penaltySVD = "ann",
                  modstr = list(gamma = gamma, nlambda = 1, lambda = lambda.seq[j]))$rank
    T.iter <- matrix(rowMeans(t(Z) - t(X) %*% B.iter), nrow = 1)
  }
  sumsq[j] <- sum((B-t(B.iter))^2)
  rank[j] <- B.rank
}




