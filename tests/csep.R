
library("free1way.docreg")

### complete separation
N <- 5
w <- gl(2, N)
y <- runif(N)
y <- c(y, y + 2)

### exact distribution
(wt <- wilcox.test(y ~ w, exact = TRUE))
(ft <- free1way(y ~ w, exact = TRUE))
stopifnot(all.equal(summary(ft, test = "Perm")$p.value,
                    wt$p.value, check.attributes = FALSE))

coef(ft)
logLik(ft)
vcov(ft)
### nonsense
confint(ft, test = "Wald")
summary(ft)
### fail with warnings
confint(ft, test = "Perm")
confint(ft, test = "Rao")
confint(ft, test = "LRT")

### asymptotic distribution
wt <- wilcox.test(y ~ w, exact = FALSE, correct = FALSE)
(ft <- free1way(y ~ w, exact = FALSE))
stopifnot(all.equal(summary(ft, test = "Perm")$p.value,
                    wt$p.value, check.attributes = FALSE))
coef(ft)
logLik(ft)
vcov(ft)
### nonsense
confint(ft, test = "Wald")
summary(ft)
### fail with warnings
confint(ft, test = "Perm")
confint(ft, test = "Rao")
confint(ft, test = "LRT")


