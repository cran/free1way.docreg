
library("free1way.docreg")

N <- 10
nd <- data.frame(y = rnorm(N), x = gl(1, N))
try(free1way(y ~ x, data = nd))

nd <- data.frame(y = rnorm(2 * N), x = gl(2, N))[seq_len(N),]
try(free1way(y ~ x, data = nd))

nd <- data.frame(y = rnorm(3 * N), x = gl(3, N))[seq_len(2 * N),]
ft1 <- free1way(y ~ x, data = nd)
ft2 <- free1way(y ~ x, data = nd)
stopifnot(all.equal(ft1, ft2))

nd <- data.frame(y = rnorm(3 * N), x = gl(3, N))
w <- rpois(nrow(nd), lambda = 2)
ft1 <- free1way(y ~ x, data = nd, weights = w)
ft2 <- free1way(y ~ x, data = nd[rep(1:nrow(nd), times = w),])
s1 <- summary(ft1)
s2 <- summary(ft2)
s1$call <- s2$call <- NULL
stopifnot(all.equal(s1, s2))

nd <- expand.grid(y = rnorm(2 * N), x = gl(2, N), block = gl(1, 2 * N))
ft1 <- free1way(y ~ x | block, data = nd)
ft2 <- free1way(y ~ x, data = nd)
s1 <- summary(ft1)
s2 <- summary(ft2)
s1$call <- s2$call <- NULL
stopifnot(all.equal(s1, s2))

