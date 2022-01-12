library(readr)
ls2_hvar_0.4b_iterB1 <- read_csv("ls2_hvar_0.4b_indepE/iterB1.csv")
ls2_hvar_0.7b_iterB1 <- read_csv("ls2_hvar_0.7b_indepE/iterB1.csv")
ls2_lvar_0.4b_iterB1 <- read_csv("ls2_lvar_0.4b_indepE/iterB1.csv")
ls2_lvar_0.7b_iterB1 <- read_csv("ls2_lvar_0.7b_indepE/iterB1.csv")

ls2_hvar_0.4b_estT <- read_csv("ls2_hvar_0.4b_indepE/estT.csv")
ls2_hvar_0.7b_estT <- read_csv("ls2_hvar_0.7b_indepE/estT.csv")
ls2_lvar_0.4b_estT <- read_csv("ls2_lvar_0.4b_indepE/estT.csv")
ls2_lvar_0.7b_estT <- read_csv("ls2_lvar_0.7b_indepE/estT.csv")

rls2_hvar_0.4b_iterB1 <- read_csv("rls2_hvar_0.4b_indepE/iterB1.csv")
rls2_hvar_0.7b_iterB1 <- read_csv("rls2_hvar_0.7b_indepE/iterB1.csv")
rls2_lvar_0.4b_iterB1 <- read_csv("rls2_lvar_0.4b_indepE/iterB1.csv")

rls3_hvar_0.4b_iterB1 <- read_csv("rls3_hvar_0.4b_indepE/iterB1.csv")
rls3_lvar_0.4b_iterB1 <- read_csv("rls3_lvar_0.4b_indepE/iterB1.csv")
rls3_hvar_0.7b_iterB1 <- read_csv("rls3_hvar_0.7b_indepE/iterB1.csv")

rrr_lvar_0.4b_iterB1 <- read_csv("rrr_lvar_0.4b_indepE/iterB1.csv")
TrueB <- read_csv("rrr_lvar_0.4b_indepE/TrueB.csv")
plot(as.matrix(rrr_lvar_0.4b_iterB1), as.matrix(TrueB))

plot(as.matrix(ls2_hvar_0.4b_estT)[1,], as.matrix(ls2_hvar_0.4b_TrueT ))
abline(a = 0, b = 1)

plot(as.matrix(ls2_lvar_0.4b_estT)[1,], as.matrix(ls2_lvar_0.4b_TrueT))
abline(a = 0, b = 1)

plot(as.matrix(ls2_hvar_0.7b_estT)[1,], as.matrix(ls2_hvar_0.7b_TrueT))
abline(a = 0, b = 1)

plot(as.matrix(ls2_lvar_0.7b_estT)[1,], as.matrix(ls2_lvar_0.7b_TrueT))
abline(a = 0, b = 1)

TrueB <- read_csv("ls2_hvar_0.4b_indepE/TrueB.csv")
ls2_hvar_0.4b_TrueT <- read_csv("ls2_hvar_0.4b_indepE/TrueT.csv")
ls2_lvar_0.4b_TrueT <- read_csv("ls2_lvar_0.4b_indepE/TrueT.csv")
ls2_hvar_0.7b_TrueT <- read_csv("ls2_hvar_0.7b_indepE/TrueT.csv")
ls2_lvar_0.7b_TrueT <- read_csv("ls2_lvar_0.7b_indepE/TrueT.csv")
X <- read_csv("ls2_hvar_0.4b_indepE/TrueB.csv")

# Changing signal to noice ratio under high variation of T(ls2 method about B)
xlabel <- expression(Estimated ~ B ~"*")
ylabel <- expression(True ~ B ~"*")
plot(as.matrix(ls2_hvar_0.4b_iterB1)[-1,], as.matrix(TrueB)[-1,],
     xlab = xlabel, ylab = ylabel, 
     cex = 0.8,
     cex.lab = 1.5,
     cex.axis = 1.5
     )
title("The performance of estimating B using 2-steps LS", adj = 0, line = 2)
mtext("under high variation of T, low signal strength", adj = 0, line = 0.5, cex = 0.8)
abline(a = 0, b = 1, lty = 2, col = 12, lwd = 2)
text(5, -4, expression(rho == 0.935), cex = 2)
cor(as.vector(as.matrix(ls2_hvar_0.4b_iterB1)[-1,]), as.vector(as.matrix(TrueB)[-1,]))

plot(as.matrix(ls2_hvar_0.7b_iterB1)[-1,], as.matrix(TrueB)[-1,],
     xlab = xlabel, ylab = ylabel, 
     cex = 0.8,
     cex.lab = 1.5,
     cex.axis = 1.5
)
title("The performance of estimating B using 2-steps LS", adj = 0, line = 2)
mtext("under high variation of T, high signal strength", adj = 0, line = 0.5, cex = 0.8)
abline(a = 0, b = 1, lty = 2, col = 12, lwd = 2)
# text(5, -4, expression(rho == 0.961), cex = 2)
text(6, -2, "Estimtion Error = 0.41", cex = 1.5)
text(6, -3, "Prediction Error = 202.32", cex = 1.5)
cor(as.vector(as.matrix(ls2_hvar_0.7b_iterB1)[-1,]), as.vector(as.matrix(TrueB)[-1,]))

plot(as.matrix(ls2_hvar_0.7b_iterB1), as.matrix(TrueB))
plot(as.matrix(ls2_hvar_0.4b_iterB1), as.matrix(ls2_hvar_0.7b_iterB1))
abline(a = 0, b = 1)

# Changing signal to noice ratio under high variation of T(ls2 method about XB)
plot(as.matrix(X) %*% t(as.matrix(ls2_hvar_0.4b_iterB1)), as.matrix(X) %*% t(as.matrix(TrueB)))
plot(as.matrix(X) %*% t(as.matrix(ls2_hvar_0.7b_iterB1)), as.matrix(X) %*% t(as.matrix(TrueB)))
plot(as.matrix(X) %*% t(as.matrix(ls2_hvar_0.4b_iterB1)), as.matrix(X) %*% t(as.matrix(ls2_hvar_0.7b_iterB1)))
abline(a = 0, b = 1)

# Changing signal to noice ratio under low variation of T(ls2 method about B)
plot(as.matrix(ls2_lvar_0.4b_iterB1), as.matrix(TrueB))
plot(as.matrix(ls2_lvar_0.7b_iterB1), as.matrix(TrueB))
plot(as.matrix(ls2_lvar_0.4b_iterB1), as.matrix(ls2_lvar_0.7b_iterB1))
abline(a = 0, b = 1)

# Changing signal to noice ratio under low variation of T(ls2 method about XB)
plot(as.matrix(X) %*% t(as.matrix(ls2_lvar_0.4b_iterB1)), as.matrix(X) %*% t(as.matrix(TrueB)))
plot(as.matrix(X) %*% t(as.matrix(ls2_lvar_0.7b_iterB1)), as.matrix(X) %*% t(as.matrix(TrueB)))
plot(as.matrix(X) %*% t(as.matrix(ls2_lvar_0.4b_iterB1)), as.matrix(X) %*% t(as.matrix(ls2_lvar_0.7b_iterB1)))
abline(a = 0, b = 1)

# Changing the variation of T under low signal to noise ratio(ls2 method about B)
plot(as.matrix(ls2_hvar_0.4b_iterB1), as.matrix(TrueB))
plot(as.matrix(ls2_lvar_0.4b_iterB1), as.matrix(TrueB))
plot(as.matrix(ls2_hvar_0.4b_iterB1), as.matrix(ls2_lvar_0.4b_iterB1))
abline(a = 0, b = 1)

# Changing the variation of T under low signal to noise ratio(ls2 method about XB)
plot(as.matrix(X) %*% t(as.matrix(ls2_hvar_0.4b_iterB1)), as.matrix(X) %*% t(as.matrix(TrueB)))
plot(as.matrix(X) %*% t(as.matrix(ls2_lvar_0.4b_iterB1)), as.matrix(X) %*% t(as.matrix(TrueB)))
plot(as.matrix(X) %*% t(as.matrix(ls2_hvar_0.4b_iterB1)), as.matrix(X) %*% t(as.matrix(ls2_lvar_0.4b_iterB1)))
abline(a = 0, b = 1)

# Changing the variation of T under high signal to noise ratio(ls2 method about B)
plot(as.matrix(ls2_hvar_0.7b_iterB1), as.matrix(TrueB))
plot(as.matrix(ls2_lvar_0.7b_iterB1), as.matrix(TrueB))
plot(as.matrix(ls2_hvar_0.7b_iterB1), as.matrix(ls2_lvar_0.7b_iterB1))

# Changing the variation of T under high signal to noise ratio(ls2 method about XB)
plot(as.matrix(X) %*% t(as.matrix(ls2_hvar_0.7b_iterB1)), as.matrix(X) %*% t(as.matrix(TrueB)))
plot(as.matrix(X) %*% t(as.matrix(ls2_lvar_0.7b_iterB1)), as.matrix(X) %*% t(as.matrix(TrueB)))
plot(as.matrix(X) %*% t(as.matrix(ls2_hvar_0.7b_iterB1)), as.matrix(X) %*% t(as.matrix(ls2_lvar_0.7b_iterB1)))
abline(a = 0, b = 1)

# The prediction error from cross validation (ls2)
ls2_hvar_0.4b_prederr <- as.matrix(read_csv("ls2_hvar_0.4b_indepE/pred_err.csv"))
ls2_hvar_0.4b_prederr_summary <- apply(ls2_hvar_0.4b_prederr, 1, min)

ls2_hvar_0.7b_prederr <- as.matrix(read_csv("ls2_hvar_0.7b_indepE/pred_err.csv"))
ls2_hvar_0.7b_prederr_summary <- apply(ls2_hvar_0.7b_prederr, 1, min)

ls2_lvar_0.4b_prederr <- as.matrix(read_csv("ls2_lvar_0.4b_indepE/pred_err.csv"))
ls2_lvar_0.4b_prederr_summary <- apply(ls2_lvar_0.4b_prederr, 1, min)

ls2_lvar_0.7b_prederr <- as.matrix(read_csv("ls2_lvar_0.7b_indepE/pred_err.csv"))
ls2_lvar_0.7b_prederr_summary <- apply(ls2_lvar_0.7b_prederr, 1, min)

boxplot(ls2_hvar_0.4b_prederr_summary,
        ls2_hvar_0.7b_prederr_summary,
        ls2_lvar_0.4b_prederr_summary,
        ls2_lvar_0.7b_prederr_summary)

# The rank estimation (ls2)
ls2_hvar_0.4b_rankest <- as.vector(as.matrix(read_csv("ls2_hvar_0.4b_indepE/rank_est.csv")))

ls2_hvar_0.7b_rankest <- as.vector(as.matrix(read_csv("ls2_hvar_0.7b_indepE/rank_est.csv")))

ls2_lvar_0.4b_rankest <- as.vector(as.matrix(read_csv("ls2_lvar_0.4b_indepE/rank_est.csv")))

ls2_lvar_0.7b_rankest <- as.vector(as.matrix(read_csv("ls2_lvar_0.7b_indepE/rank_est.csv")))


rls2_hvar_0.4b_rankest <- as.vector(as.matrix(read_csv("rls2_hvar_0.4b_indepE/rank_est.csv")))

rls2_hvar_0.7b_rankest <- as.vector(as.matrix(read_csv("rls2_hvar_0.7b_indepE/rank_est.csv")))


rls3_hvar_0.4b_rankest <- as.vector(as.matrix(read_csv("rls3_hvar_0.4b_indepE/rank_est.csv")))

rls3_hvar_0.7b_rankest <- as.vector(as.matrix(read_csv("rls3_hvar_0.7b_indepE/rank_est.csv")))


rrr_hvar_0.4b_rankest <- as.vector(as.matrix(read_csv("rrr_hvar_0.4b_indepE/rank_est.csv")))

rrr_hvar_0.7b_rankest <- as.vector(as.matrix(read_csv("rrr_hvar_0.7b_indepE/rank_est.csv")))


boxplot(ls2_hvar_0.4b_rankest,
        ls2_hvar_0.7b_rankest,
        rls2_hvar_0.4b_rankest,
        rls2_hvar_0.7b_rankest, 
        rls3_hvar_0.4b_rankest,
        rls3_hvar_0.7b_rankest,
        rrr_hvar_0.4b_rankest,
        rrr_hvar_0.7b_rankest,
        names = c("LS2, lb", "LS2, hb", "RLS2, lb", "RLS2, hb", "RLS3, lb", "RLS3, hb", "RRR, lb", "RRR, hb"))

# The rank estimation (ls2)
ls2_hvar_0.4b_congiter <- as.vector(as.matrix(read_csv("ls2_hvar_0.4b_indepE/cong_iter.csv")))

ls2_hvar_0.7b_congiter <- as.vector(as.matrix(read_csv("ls2_hvar_0.7b_indepE/cong_iter.csv")))


rls2_hvar_0.4b_congiter <- as.vector(as.matrix(read_csv("rls2_hvar_0.4b_indepE/cong_iter.csv")))

rls2_hvar_0.7b_congiter <- as.vector(as.matrix(read_csv("rls2_hvar_0.7b_indepE/cong_iter.csv")))


rls3_hvar_0.4b_congiter <- as.vector(as.matrix(read_csv("rls3_hvar_0.4b_indepE/cong_iter.csv")))

rls3_hvar_0.7b_congiter <- as.vector(as.matrix(read_csv("rls3_hvar_0.7b_indepE/cong_iter.csv")))




boxplot(ls2_hvar_0.4b_congiter,
        ls2_hvar_0.7b_congiter,
        rls2_hvar_0.4b_congiter,
        rls2_hvar_0.7b_congiter, 
        rls3_hvar_0.4b_congiter,
        rls3_hvar_0.7b_congiter,
        names = c("LS2, lb", "LS2, hb", "RLS2, lb", "RLS2, hb", "RLS3, lb", "RLS3, hb"))

# Changing the variation of T (rls2 method about B)

plot(as.matrix(rls2_hvar_0.4b_iterB1)[-1,], as.matrix(TrueB)[-1,],
     xlab = xlabel, ylab = ylabel, 
     cex = 0.8,
     cex.lab = 1.5,
     cex.axis = 1.5
)
title("The performance of estimating B using 2-steps RLS", adj = 0, line = 2)
mtext("under high variation of T, low signal strength", adj = 0, line = 0.5, cex = 0.8)
abline(a = 0, b = 1, lty = 2, col = 12, lwd = 2)
text(4, -4, expression(rho == 0.877), cex = 2)
cor(as.vector(as.matrix(rls2_hvar_0.4b_iterB1)[-1,]), as.vector(as.matrix(TrueB)[-1,]))

plot(as.matrix(rls2_hvar_0.7b_iterB1)[-1,], as.matrix(TrueB)[-1,],
     xlab = xlabel, ylab = ylabel, 
     cex = 0.8,
     cex.lab = 1.5,
     cex.axis = 1.5
)
title("The performance of estimating B using 2-steps RLS", adj = 0, line = 2)
mtext("under high variation of T, high signal strength", adj = 0, line = 0.5, cex = 0.8)
abline(a = 0, b = 1, lty = 2, col = 12, lwd = 2)
text(7, -4, expression(rho == 0.902), cex = 2)
cor(as.vector(as.matrix(rls2_lvar_0.4b_iterB1)[-1,]), as.vector(as.matrix(TrueB)[-1,]))

plot(as.matrix(rls2_hvar_0.4b_iterB1), as.matrix(TrueB))
plot(as.matrix(rls2_lvar_0.4b_iterB1), as.matrix(TrueB))
plot(as.matrix(rls2_hvar_0.4b_iterB1), as.matrix(rls2_lvar_0.4b_iterB1))

# Changing the variation of T (rls2 method about XB)
plot(as.matrix(X) %*% t(as.matrix(rls2_hvar_0.4b_iterB1)), as.matrix(X) %*% t(as.matrix(TrueB)))
plot(as.matrix(X) %*% t(as.matrix(rls2_lvar_0.4b_iterB1)), as.matrix(X) %*% t(as.matrix(TrueB)))
plot(as.matrix(X) %*% t(as.matrix(rls2_hvar_0.4b_iterB1)), as.matrix(X) %*% t(as.matrix(rls2_lvar_0.4b_iterB1)))

# The prediction error from cross validation (rls2)
rls2_hvar_0.4b_prederr <- as.matrix(read_csv("rls2_hvar_0.4b_indepE/pred_err.csv"))
rls2_hvar_0.4b_prederr_summary <- apply(rls2_hvar_0.4b_prederr, 1, min)

rls2_hvar_0.7b_prederr <- as.matrix(read_csv("rls2_hvar_0.7b_indepE/pred_err.csv"))
rls2_hvar_0.7b_prederr_summary <- apply(rls2_hvar_0.7b_prederr, 1, min)

rls2_lvar_0.4b_prederr <- as.matrix(read_csv("rls2_lvar_0.4b_indepE/pred_err.csv"))
rls2_lvar_0.4b_prederr_summary <- apply(rls2_lvar_0.4b_prederr, 1, min)

rls2_lvar_0.7b_prederr <- as.matrix(read_csv("rls2_lvar_0.7b_indepE/pred_err.csv"))
rls2_lvar_0.7b_prederr_summary <- apply(rls2_lvar_0.7b_prederr, 1, min)

boxplot(rls2_hvar_0.4b_prederr_summary,
        rls2_lvar_0.4b_prederr_summary)

# Rank estimation (rls2)
rls2_hvar_0.4b_rankest <- as.matrix(read_csv("rls2_hvar_0.4b_indepE/rank_est.csv"))
rls2_lvar_0.4b_rankest <- as.matrix(read_csv("rls2_lvar_0.4b_indepE/rank_est.csv"))
mean(rls2_hvar_0.4b_rankest)
mean(rls2_lvar_0.4b_rankest)

# Changing the variation of T (rls3 method about B)
plot(as.matrix(rls3_hvar_0.4b_iterB1), as.matrix(TrueB)[-1,])
plot(as.matrix(rls3_lvar_0.4b_iterB1), as.matrix(TrueB)[-1,])
plot(as.matrix(rls3_hvar_0.4b_iterB1), as.matrix(rls3_lvar_0.4b_iterB1))

# Changing the variation of T (rls3 method about XB)
plot(as.matrix(X) %*% t(as.matrix(rls3_hvar_0.4b_iterB1)), as.matrix(X) %*% t(as.matrix(TrueB)[-1,]))
plot(as.matrix(X) %*% t(as.matrix(rls3_lvar_0.4b_iterB1)), as.matrix(X) %*% t(as.matrix(TrueB)[-1,]))
plot(as.matrix(X) %*% t(as.matrix(rls3_hvar_0.4b_iterB1)), as.matrix(X) %*% t(as.matrix(rls3_lvar_0.4b_iterB1)))

# Changing the signal to noise ratio (rls3 method about B)

plot(as.matrix(rls3_hvar_0.4b_iterB1), as.matrix(TrueB)[-1,],
     xlab = xlabel, ylab = ylabel, 
     cex = 0.8,
     cex.lab = 1.5,
     cex.axis = 1.5
)
title("The performance of estimating B using 3-steps RLS", adj = 0, line = 2)
mtext("under high variation of T, low signal strength", adj = 0, line = 0.5, cex = 0.8)
abline(a = 0, b = 1, lty = 2, col = 12, lwd = 2)
text(3, -4, expression(rho == 0.870), cex = 2)
cor(as.vector(as.matrix(rls3_hvar_0.4b_iterB1)), as.vector(as.matrix(TrueB)[-1,]))

plot(as.matrix(rls3_hvar_0.7b_iterB1), as.matrix(TrueB)[-1,],
     xlab = xlabel, ylab = ylabel, 
     cex = 0.8,
     cex.lab = 1.5,
     cex.axis = 1.5
)
title("The performance of estimating B using 3-steps RLS", adj = 0, line = 2)
mtext("under high variation of T, high signal strength", adj = 0, line = 0.5, cex = 0.8)
abline(a = 0, b = 1, lty = 2, col = 12, lwd = 2)
text(6, -4, expression(rho == 0.889), cex = 2)
cor(as.vector(as.matrix(rls3_hvar_0.7b_iterB1)), as.vector(as.matrix(TrueB)[-1,]))

plot(as.matrix(rls3_hvar_0.4b_iterB1), as.matrix(TrueB)[-1,])
plot(as.matrix(rls3_hvar_0.7b_iterB1), as.matrix(TrueB)[-1,])
plot(as.matrix(rls3_hvar_0.4b_iterB1), as.matrix(rls3_hvar_0.7b_iterB1))

# Changing the signal to noise ratio (rls3 method about XB)
plot(as.matrix(rls3_hvar_0.4b_iterB1), as.matrix(TrueB)[-1,])
plot(as.matrix(rls3_hvar_0.7b_iterB1), as.matrix(TrueB)[-1,])
plot(as.matrix(rls3_hvar_0.4b_iterB1), as.matrix(rls3_hvar_0.7b_iterB1))

# The prediction error from cross validation (rls3)
rls3_hvar_0.4b_prederr <- as.matrix(read_csv("rls3_hvar_0.4b_indepE/pred_err.csv"))
rls3_hvar_0.4b_prederr_summary <- apply(rls3_hvar_0.4b_prederr, 1, min)

rls3_hvar_0.7b_prederr <- as.matrix(read_csv("rls3_hvar_0.7b_indepE/pred_err.csv"))
rls3_hvar_0.7b_prederr_summary <- apply(rls3_hvar_0.7b_prederr, 1, min)

rls3_lvar_0.4b_prederr <- as.matrix(read_csv("rls3_lvar_0.4b_indepE/pred_err.csv"))
rls3_lvar_0.4b_prederr_summary <- apply(rls3_lvar_0.4b_prederr, 1, min)

rls3_lvar_0.7b_prederr <- as.matrix(read_csv("rls3_lvar_0.7b_indepE/pred_err.csv"))
rls3_lvar_0.7b_prederr_summary <- apply(rls3_lvar_0.7b_prederr, 1, min)

boxplot(rls3_hvar_0.4b_prederr_summary,
        rls3_lvar_0.4b_prederr_summary)

# Rank estimation (rls3)
rls3_hvar_0.4b_rankest <- as.matrix(read_csv("rls3_hvar_0.4b_indepE/rank_est.csv"))
rls3_lvar_0.4b_rankest <- as.matrix(read_csv("rls3_lvar_0.4b_indepE/rank_est.csv"))
mean(rls3_hvar_0.4b_rankest)
mean(rls3_lvar_0.4b_rankest)

# The prediction error from cross validation (rrr)
rrr_hvar_0.4b_prederr <- as.matrix(read_csv("rrr_hvar_0.4b_indepE/pred_err.csv"))
rrr_hvar_0.4b_prederr_summary <- apply(rrr_hvar_0.4b_prederr, 1, min)

rrr_hvar_0.7b_prederr <- as.matrix(read_csv("rrr_hvar_0.7b_indepE/pred_err.csv"))
rrr_hvar_0.7b_prederr_summary <- apply(rrr_hvar_0.7b_prederr, 1, min)

rrr_lvar_0.4b_prederr <- as.matrix(read_csv("rrr_lvar_0.4b_indepE/pred_err.csv"))
rrr_lvar_0.4b_prederr_summary <- apply(rrr_lvar_0.4b_prederr, 1, min)

rrr_lvar_0.7b_prederr <- as.matrix(read_csv("rrr_lvar_0.7b_indepE/pred_err.csv"))
rrr_lvar_0.7b_prederr_summary <- apply(rrr_lvar_0.7b_prederr, 1, min)

boxplot(rrr_hvar_0.4b_prederr_summary,
        rrr_hvar_0.7b_prederr_summary,
        rrr_lvar_0.4b_prederr_summary,
        rrr_lvar_0.7b_prederr_summary)

rrr_hvar_0.4b_iterB1 <- read_csv("rrr_hvar_0.4b_indepE/iterB1.csv")
rrr_hvar_0.7b_iterB1 <- read_csv("rrr_hvar_0.7b_indepE/iterB1.csv")
rrr_lvar_0.4b_iterB1 <- read_csv("rrr_lvar_0.4b_indepE/iterB1.csv")
rrr_lvar_0.7b_iterB1 <- read_csv("rrr_lvar_0.7b_indepE/iterB1.csv")
plot(as.matrix(rrr_hvar_0.4b_iterB1)[-1,], as.matrix(TrueB)[-1,],
     xlab = xlabel, ylab = ylabel, 
     cex = 0.8, 
     cex.lab = 1.5,
     cex.axis = 1.5
)
title("The performance of estimating B using RRR", adj = 0, line = 2.5, cex.main = 2.5)
# title("The performance of estimating B using RRR")
mtext("under high variation of T, low signal strength", adj = 0, line = 0.5, cex = 2)
abline(a = 0, b = 1, lty = 2, col = 12, lwd = 2)
text(-20, 5, expression(rho == -0.011), cex = 2)
cor(as.vector(as.matrix(rrr_hvar_0.4b_iterB1)[-1,]), as.vector(as.matrix(TrueB)[-1,]))

plot(as.matrix(rrr_hvar_0.7b_iterB1)[-1,], as.matrix(TrueB)[-1,],
     xlab = xlabel, ylab = ylabel, 
     cex = 0.8,
     cex.lab = 1.5,
     cex.axis = 1.5
)
title("The performance of estimating B using RRR", adj = 0, line = 2)
mtext("under high variation of T, high signal strength", adj = 0, line = 0.5, cex = 0.8)
abline(a = 0, b = 1, lty = 2, col = 12, lwd = 2)
text(-20, 5, expression(rho == 0.063), cex = 2)
cor(as.vector(as.matrix(rrr_hvar_0.7b_iterB1)[-1,]), as.vector(as.matrix(TrueB)[-1,]))


attach(mtcars)
par(mfrow=c(1,2))
plot(as.matrix(ls2_hvar_0.7b_iterB1)[-1,], as.matrix(TrueB)[-1,],
     xlab = xlabel, ylab = ylabel, 
     cex = 0.8,
     cex.lab = 1.5,
     cex.axis = 1.5
)
title("2-steps LS", adj = 0, line = 2, cex = 1)
mtext("under high variation of T, high signal strength \n Estimation Error = 0.41, Prediction Error = 202.32", cex = 0.7)
abline(a = 0, b = 1, lty = 2, col = 12, lwd = 2)
# text(5, -4, expression(rho == 0.961), cex = 2)

plot(as.matrix(rrr_hvar_0.7b_iterB1)[-1,], as.matrix(TrueB)[-1,],
     xlab = xlabel, ylab = ylabel, 
     cex = 0.8,
     cex.lab = 1.5,
     cex.axis = 1.5
)
title("RRR", adj = 0, line = 2, cex = 1)
mtext("under high variation of T, high signal strength \n Estimation Error = 90.04, Prediction Error = 2293.71", cex = 0.7)
abline(a = 0, b = 1, lty = 2, col = 12, lwd = 2)
# text(-20, 5, expression(rho == 0.063), cex = 2)
