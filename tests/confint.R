
library("free1way.docreg")
options(digits = 5)
set.seed(29)
tol <- 1e-5

level <- .9
N <- 15
w <- gl(2, N)
beta <- -1
y <- rlogis(length(w), location = c(0, beta)[w])

ft <- free1way(y ~ w)
coef(ft)
logLik(ft)
(ci <- confint(ft, test = "Wald", level = level))
summary(ft1 <- free1way(y ~ w, mu = ci[1]), test = "Wald", alternative =
"greater")$p.value
summary(ft2 <- free1way(y ~ w, mu = ci[2]), test = "Wald", alternative =
"less")$p.value
summary(ft12 <- free1way(y ~ w, mu = ci[1]), test = "Wald")$p.value
summary(ft22 <- free1way(y ~ w, mu = ci[2]), test = "Wald")$p.value
coef(ft1) + ci[1]
coef(ft2) + ci[2]
logLik(ft1)
logLik(ft2)

all.equal(ci, confint(ft1, test = "Wald", level = level) + ci[1], 
          tolerance = tol)
all.equal(ci, confint(ft2, test = "Wald", level = level) + ci[2],
          tolerance = tol)
all.equal(ci, confint(ft12, test = "Wald", level = level) + ci[1],
          tolerance = tol)
all.equal(ci, confint(ft22, test = "Wald", level = level) + ci[2],
          tolerance = tol)


(ci <- confint(ft, test = "LRT", level = level))
summary(ft12 <- free1way(y ~ w, mu = ci[1]), test = "LRT")$p.value
summary(ft22 <- free1way(y ~ w, mu = ci[2]), test = "LRT")$p.value

all.equal(ci, confint(ft12, test = "LRT", level = level) + ci[1],
          tolerance = tol)
all.equal(ci, confint(ft22, test = "LRT", level = level) + ci[2],
          tolerance = tol)


(ci <- confint(ft, test = "Rao", level = level))
summary(ft12 <- free1way(y ~ w, mu = ci[1]), test = "Rao")$p.value
summary(ft22 <- free1way(y ~ w, mu = ci[2]), test = "Rao")$p.value

all.equal(ci, confint(ft12, test = "Rao", level = level) + ci[1],
          tolerance = tol)
all.equal(ci, confint(ft22, test = "Rao", level = level) + ci[2],
          tolerance = tol)


(ci <- confint(ft, test = "Permutation", level = level))
summary(ft1 <- free1way(y ~ w, mu = ci[1]), test = "Permutation", alternative = "greater")$p.value
summary(ft2 <- free1way(y ~ w, mu = ci[2]), test = "Permutation", alternative = "less")$p.value
summary(ft12 <- free1way(y ~ w, mu = ci[1]), test = "Permutation")$p.value
summary(ft22 <- free1way(y ~ w, mu = ci[2]), test = "Permutation")$p.value
coef(ft1) + ci[1]
coef(ft2) + ci[2]
logLik(ft1)
logLik(ft2)

all.equal(ci, confint(ft1, test = "Permutation", level = level) + ci[1], 
          tolerance = tol)
all.equal(ci, confint(ft2, test = "Permutation", level = level) + ci[2],
          tolerance = tol)
all.equal(ci, confint(ft12, test = "Permutation", level = level) + ci[1],
          tolerance = tol)
all.equal(ci, confint(ft22, test = "Permutation", level = level) + ci[2],
          tolerance = tol)

ftc <- free1way(y ~ w, mu = coef(ft))
summary(ftc, test = "Permutation", alternative = "greater")$p.value
summary(ftc, test = "Permutation", alternative = "less")$p.value
summary(ftc, test = "Permutation")$p.value

confint(ftc, test = "Permutation", level = level) + coef(ft)

B <- 1000
set.seed(29)
ft <- free1way(y ~ w, B = B)
(ci <- confint(ft, test = "Permutation", level = level))
### check inversion of confidence intervals: These numbers should be
set.seed(29)
### <= (1 - level) / 2 = .05
summary(ft1 <- free1way(y ~ w, B = B, mu = ci[1]), test = "Permutation", alternative = "greater")$p.value
set.seed(29)
summary(ft2 <- free1way(y ~ w, B = B, mu = ci[2]), test = "Permutation", alternative = "less")$p.value
set.seed(29)
### <= 1 - level = .1
summary(ft12 <- free1way(y ~ w, B = B, mu = ci[1]), test = "Permutation")$p.value
set.seed(29)
summary(ft22 <- free1way(y ~ w, B = B, mu = ci[2]), test = "Permutation")$p.value
coef(ft1) + ci[1]
coef(ft2) + ci[2]
logLik(ft1)
logLik(ft2)

all.equal(ci, confint(ft1, test = "Permutation", level = level) + ci[1], 
          tolerance = tol)
all.equal(ci, confint(ft2, test = "Permutation", level = level) + ci[2],
          tolerance = tol)
all.equal(ci, confint(ft12, test = "Permutation", level = level) + ci[1],
          tolerance = tol)
all.equal(ci, confint(ft22, test = "Permutation", level = level) + ci[2],
          tolerance = tol)


summary(ftg <- free1way(y ~ w, mu = coef(ft)), test = "Permutation", alternative = "greater")$p.value
summary(ftl <- free1way(y ~ w, mu = coef(ft)), test = "Permutation", alternative = "less")$p.value
summary(ftt <- free1way(y ~ w, mu = coef(ft)), test = "Permutation")$p.value

set.seed(29)
ftc <- free1way(y ~ w, mu = coef(ft), B = B)
summary(ftc, test = "Permutation", alternative = "greater")$p.value
summary(ftc, test = "Permutation", alternative = "less")$p.value
summary(ftc, test = "Permutation")$p.value

confint(ftc, test = "Permutation", level = level) + coef(ft)
