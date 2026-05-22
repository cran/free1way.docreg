
library("free1way.docreg")
options(digits = 4)

args <- expand.grid(n = c(15, 50),
                    delta = c(0, 1),
                    alloc_ratio = c(1, 3))
args$seed <- 290875
args$pwr <- 0
p <- which(names(args) == "pwr")

for (i in 1:nrow(args))
    args$pwr[i] <- do.call("power.free1way.test", args[i,-p])$power

args

### SWIFT trial
delta <- log(1.38)
x <- power.free1way.test(delta = delta, power = .8, seed = 29)
x1 <- power.free1way.test(n = x$n, power = x$power, delta = delta, sig.level = NULL, seed = 29)
x2 <- power.free1way.test(n = x$n, delta = delta, seed = 29)
x3 <- power.free1way.test(n = x$n, power = x$power, seed = 29)

stopifnot(isTRUE(all.equal(x, x1, tol = 1e-2)))
stopifnot(isTRUE(all.equal(x, x2, tol = 1e-2)))
stopifnot(isTRUE(all.equal(x, x3, tol = 1e-2)))


