
pkgs <- c("free1way.docreg", "rms")
if (!all(sapply(pkgs, require, character.only = TRUE)))
    quit()

set.seed(2908)

### speed comparisons for larger sample sizes
N <- 20000
K <- floor(N / 10)

### categorical, rfree1way is faster
d <- rfree1way(n = N, delta = .25, prob = rep(1 / K, K))
system.time(m0 <- free1way(y ~ groups, data = d))
system.time(m1 <- orm(y ~ groups, data = d))

stopifnot(isTRUE(all.equal(c(logLik(m0)), c(logLik(m1)))))

stopifnot(isTRUE(all.equal(coef(m0), rev(coef(m1))[1], check.attributes = FALSE)))

stopifnot(isTRUE(all.equal(c(vcov(m0)), vcov(m1)[2,2], check.attributes = FALSE)))

### numerical, orm is faster
d <- rfree1way(n = N, delta = .25)
d$y <- recode2integer(d$y)$y ### rounding as in orm internally
system.time(m0 <- free1way(y ~ groups, data = d))
system.time(m1 <- orm(y ~ groups, data = d))

stopifnot(isTRUE(all.equal(c(logLik(m0)), c(logLik(m1)))))

stopifnot(isTRUE(all.equal(coef(m0), rev(coef(m1))[1], check.attributes = FALSE)))

stopifnot(isTRUE(all.equal(c(vcov(m0)), vcov(m1)[2,2], check.attributes = FALSE)))

### effect of rounding
length(coef(m1))
