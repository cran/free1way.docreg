
pkgs <- c("free1way.docreg", "rms", "survival")
if (!all(sapply(pkgs, require, character.only = TRUE, warn.conflicts = FALSE)))
    quit()

options(digits = 5)
set.seed(29)

N <- 500
x <- gl(2, N)
y <- rlogis(length(x), location = c(0, 2)[x])

ci <- confint(m <- orm(y ~ x))
c(rev(coef(m))[1], ci[nrow(ci),])
logLik(m)

d <- data.frame(y = y, x = x)

(ft <- free1way(y = y, groups = x))
coef(ft)
confint(ft)
logLik(ft)
(ft <- free1way(y ~ x, data = d))
(ft <- free1way(y = y, groups = x, nbins = 99))
(ft <- free1way(y = y, groups = x, nbins = 100))

-coef(m <- coxph(Surv(y, rep(TRUE, length(y))) ~ x))
-rev(confint(m))
### rms: cloglog for P(Y < y) is loglog for P(Y >= y), the latter
### being used in orm()
ci <- confint(m <- orm(y ~ x, family = "loglog"))
c(rev(coef(m))[1], ci[nrow(ci),])
logLik(m)

ft <- free1way(y = y, groups = x, link = "cloglog")
coef(ft)
confint(ft)
logLik(ft)
ft <- free1way(y = y, groups = x, link = "cloglog", nbins = 99)
coef(ft)
confint(ft)
ft <- free1way(y = y, groups = x, link = "cloglog", nbins = 100)
coef(ft)
confint(ft)

ci <- confint(m <- orm(y ~ x, family = "cloglog"))
c(rev(coef(m))[1], ci[nrow(ci),])
logLik(m)

ft <- free1way(y = y, groups = x, link = "loglog")
coef(ft)
confint(ft)
logLik(ft)
ft <- free1way(y = y, groups = x, link = "loglog", nbins = 99)
coef(ft)
confint(ft)
ft <- free1way(y = y, groups = x, link = "loglog", nbins = 100)
coef(ft)
confint(ft)

### probit
ci <- confint(m <- orm(y ~ x, family = "probit"))
c(rev(coef(m))[1], ci[nrow(ci),])
logLik(m)

ft <- free1way(y = y, groups = x, link = "probit")
coef(ft)
confint(ft)
logLik(ft)
ft <- free1way(y = y, groups = x, link = "probit", nbins = 99)
coef(ft)
confint(ft)
ft <- free1way(y = y, groups = x, link = "probit", nbins = 100)
coef(ft)
confint(ft)


### with offset
mu <- 1
off <- (x == levels(x)[2]) * mu

ci <- confint(m <- orm(y ~ x + offset(off)))
c(rev(coef(m))[1], ci[nrow(ci),])
logLik(m)

ft <- free1way(y = y, groups = x, mu = 1)
coef(ft)
confint(ft)
logLik(m)
ft <- free1way(y = y, groups = x, mu = 1, nbins = 99)
coef(ft)
confint(ft)
ft <- free1way(y = y, groups = x, mu = 1, nbins = 100)
coef(ft)
confint(ft)

### permutations
N <- 15
x <- gl(2, N)
y <- rlogis(length(x), location = c(0, 2)[x])

ci <- confint(m <- orm(y ~ x))
c(rev(coef(m))[1], ci[nrow(ci),])
logLik(m)

ft <- free1way(y = y, groups = x)
print(ft, test = "Wald")
coef(ft)
confint(ft, test = "Wald")
logLik(m)
print(ft, test = "Rao")
confint(ft, test = "Rao")
print(ft, test = "LRT")
confint(ft, test = "LRT")
print(ft, test = "Permutation")
confint(ft, test = "Permutation")

ft <- free1way(y = y, groups = x, B = 10000)
print(ft, test = "Permutation")
confint(ft, test = "Permutation")

### generate data from linear transformation model with linear predictor
### as offset in rfree1way
N <- 500
K <- 3
P <- 2
beta <- c(-1.5, 1.75)
X <- matrix(runif(N * K * P), ncol = P)
lp <- X %*% beta
d <- rfree1way(n = N, delta = c(0, .5), offset = lp)
d <- cbind(d, as.data.frame(X))
tail(coef(orm(y ~ groups + V1 + V2, data = d)), 4)
