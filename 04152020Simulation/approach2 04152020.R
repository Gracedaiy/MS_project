source("data simulation.R")
library(rrpack)
library(MASS)
library(matrixcalc)
##functions
rrr <- function(X,Y,lambda,gamma = 2){
  rrresult <- list()
  #1 estimate Clhat
  Clhat <- ginv(t(X)%*%X)%*%t(X)%*% Y
  
  #2 eigendecomposite of CLhar
  Uhat <- svd(X%*%Clhat)$u
  vhat <- svd(X%*%Clhat)$v
  D <- diag(svd(X%*%Clhat)$d)
  
  #3. weight for soft-thredshold and rank
  weight = (svd(X%*%Clhat)$d)^-2
  rrresult$r = sum(svd(X%*%Clhat)$d > lambda^(1/(2+1)))
  
  #4. soft-thredsholding
  diagvalue = svd(X%*%Clhat)$d -lambda*weight
  th <- ifelse(diagvalue < 0, 0, diagvalue)
  SDhat = diag(th)
  
  rrresult$thd <- diagvalue
  
  #5. low rank solution C
  rrresult$CShat = Clhat %*% vhat %*% ginv(D) %*% SDhat %*% t(vhat)
  
  return(rrresult)
}


##initial beta0, B, and T, lambda
simplealgorithm <- function(X,W,Y, eb,ey,et,ewb,ebeta0,lambda,simN){
  result <- list()
  
  p <- dim(X)[2]
  q <- dim(Y)[2]
  beta0hat <- apply(Y, 2, mean)
  Bhat <- matrix(1, nrow = p, ncol = q)
  That <- log(apply(X, 1, sum))
  WBhat <- W %*% Bhat
  Yhat <- beta0hat +  W%*%Bhat + That %*% t(Ip)%*%Bhat
  
  #initial convergence
  convergeb = F
  convergey = F
  convergebeta0 = F
  convergewb = F
  converget = F
  
  # 1 vector
  In <- rep(1,n)
  Ip <- rep(1,p)
  iter = 1
  
  ## repeat until converge
  diffb = 10
  diffy = 10
  difft = 10
  diffWBhat = 10
  diffbeta0 = 10
  
  while(((convergeb == F || convergebeta0 == F) && (convergey == F || converget == F || convergewb == F)) && iter <= simN){
  
  #start k+1 step
  #1.fixed T, optimize g over beta0 and Bhatnew
    X = cbind(rep(1,n), W + That %*% t(Ip))
    rrresult <-  rrr(X = X, Y = Y, lambda = lambda, gamma = 1)
    CShat <- rrresult$CShat
    
  #2.fixed beta0 and B, optimize g over T
    Z = cbind(rep(1,n), W)
    S = Y - Z %*% CShat
    
    Bhatnew <- CShat[2:(p+1),]
    beta0hatnew <- CShat[1,]
    
    Thatnew = t(solve(Ip %*% Bhatnew %*% t(Bhatnew)%*%Ip) %*% t(t(Bhatnew)%*% Ip) %*% t(S))
    #A = svd(S)$u[,1] %*%t(svd(S)$v[,1])*svd(S)$d[1]
 
  #3.convergence? by Bhat and beta0 and WB
    diffb <- max(abs(Bhat - Bhatnew))
    if (diffb < eb)
    {convergeb = T}
    else if(diffb > eb){
      convergeb = F
    }
    
    diffbeta0 <- max(abs(beta0hat - beta0hatnew))
    if (diffbeta0 < ebeta0)
    {convergebeta0 = T}
    else if(diffbeta0 > ebeta0){
      convergebeta0 = F
    }
    
    diffWBhat <- max(abs(WBhat - W%*%Bhatnew))
    if (diffWBhat < ewb)
    {convergewb = T}
    else if(diffWBhat > ewb){
      convergewb = F
    }
    
  #4. convergence? by y-yhat
    Yhatnew <- beta0hatnew +  W%*%Bhatnew + That %*% t(Ip)%*%Bhatnew
    diffy <- mean(abs(Yhat-Yhatnew))
    if (diffy < ey)
    {convergey = T}
    else if(diffy > ey){
      convergey = F
      Bhat <- Bhatnew
    }
    iter = iter + 1
    print(iter)
    
  #5. Congernce? by T
    difft <- max(abs(That - Thatnew))
    if (difft < et)
    {converget = T}
    else if(difft > et){
      converget = F
    } 
  
  ##result 
    beta0hat <- beta0hatnew
    Bhat <- Bhatnew
    WBhat <-  W%*%Bhatnew
    That <- Thatnew
    Yhat<- Yhatnew
  }


  #6. difference in objective function
  result$convergeb <- convergeb
  result$convergey <- convergey
  result$convergebeta0 <- convergebeta0
  result$convergewb <-convergewb 
  result$converget <- converget
  
  result$diffy <- diffy
  result$diffb <- diffb
  result$difft <- difft
  result$diffbeta0 <- diffbeta0
  result$diffWBat <- diffWBhat
  
  result$iter <- iter
  result$r <- rrresult$r
  
  result$coef <- Bhatnew
  result$That <- That
  result$thd <- rrresult$thd
  result$beta0hat <- beta0hat
  result$Yhat <- Yhat
  
return(result)
}

###simulate data

p = 100
q = 50
n = 200
set.seed(333)
utrue <- t(replicate(p, rnorm(q, mean = 0, sd = 1) ,simplify = "matrix"))
dtrue <- diag(c(10:1, rep(0, q-10)))
vtrue <- t(replicate(q, rnorm(q, mean = 0, sd = 1) ,simplify = "matrix"))

Btrue <- utrue %*% dtrue %*% t(vtrue)
microdata <- microbsimudata(n = 200,p = 100, q = 50, B = Btrue)
str(microdata)
Y <- microdata$Y
X <- t(microdata$X)
W <- microdata$W
dim(Y)
dim(X)
dim(W)

Bsimulation <- simplealgorithm(X =X, W = W, Y = Y, eb = 10^-9, ey = 10^-9,ebeta0 = 10^-9,ewb = 10^-9,et = 10^-9, lambda = 3300, simN = 500)

Bsimulation$r
Bsimulation$diffy
Bsimulation$diffb
Bsimulation$difft
Bsimulation$diffWBat
Bsimulation$diffbeta0

Bsimulation$convergeb 
Bsimulation$convergey 
Bsimulation$convergewb 
Bsimulation$convergebeta0
Bsimulation$converget

Bsimulation$That
Bsimulation$beta0hat
View(Bsimulation$coef - Btrue)
View(Bsimulation$That - log(xt))
View(Bsimulation$Yhat - Y)
View(W %*% Bsimulation$coef - W%*%Btrue)

temp = cbind(rep(1,n), W + That %*% t(Ip))
B <- rbind(Bsimulation$beta0hat, Bsimulation$coef)
Btruewithbeta0 <- rbind(rep(1,q), Btrue)
hist(temp %*% B- temp%*% Btruewithbeta0,breaks = 50)
hist(B- Btruewithbeta0,breaks = 50)
hist(Bsimulation$coef- Btrue,breaks = 50)

svd(temp %*% B)$d
svd(temp%*% Btruewithbeta0)$d

svd(temp %*% B)$v[1:10,1:10]
svd(temp%*% Btruewithbeta0)$u[1:10,1:10]

svd(Bsimulation$coef)$d
svd(Btrue)$d

svd(Bsimulation$coef)$v[1:10,1:10]
svd(Btrue)$v[1:10,1:10]
