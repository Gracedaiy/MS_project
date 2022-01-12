X <- as.matrix(X)
B1_1 <- as.matrix(B1_1)
B6_1 <- as.matrix(B6_1)
iterB1 <- as.matrix(iterB1)
TrueB <- as.matrix(TrueB)
iterT1 <- as.matrix(iterT1)

sum(abs(B1_1 - iterB1))/sum(abs(iterB1))
sum((X %*% B1_1 + matrix(rep(iterT1, each = 40), nrow = 100, byrow = T) - 
      X %*% iterB1 - matrix(rep(5.38, each = 4000), nrow = 100, byrow = T))^2)
sum(abs(X %*% TrueB - X %*% iterB1)^2)
sum((X %*% TrueB - X %*% B1_1)^2)


png(paste0("estB1.png"))
plot(B1_1, iterB1)
abline(a = 0, b = 1)
dev.off()

png(paste0("estXB1.png"))
plot(X %*% B1_1, X %*% iterB1)
abline(a = 0, b = 1)
dev.off()

plot(iterT1)

sum((iterT1 - 5.38)^2)
sum((B1_1 - iterB1)^2)
sum((B2_1 - iterB1)^2)
sum((X %*% B1_1 - X %*% iterB1)^2)
sum((X %*% B2_1 - X %*% iterB1)^2)
sum((X %*% B2_1 + matrix(rep(mean(iterT1), each = 4000), nrow = 100, byrow = T) - 
       X %*% iterB1 - matrix(rep(5.38, each = 4000), nrow = 100, byrow = T))^2)

plot(B1_1, TrueB)
plot(iterB1, TrueB)
sum((B1_1 - TrueB)^2)
sum((iterB1 - TrueB)^2)

sum((X %*% B1_1 - X %*% TrueB)^2)
sum((X %*% iterB1 - X %*% TrueB)^2)



plot(B6_1, TrueB)
plot(iterB1, TrueB)
sum((B6_1 - TrueB)^2)
sum((iterB1 - TrueB)^2)

sum((X %*% B6_1 - X %*% TrueB)^2)
sum((X %*% iterB1 - X %*% TrueB)^2)


sum((X %*% B6_1 - X %*% iterB1)^2)
