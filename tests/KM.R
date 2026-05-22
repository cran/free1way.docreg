
### check blockwise parameterisation and null parameter estimation
### against Kaplan-Meier

pkgs <- c("free1way.docreg", "survival")
if (!all(sapply(pkgs, require, character.only = TRUE)))
    quit()

set.seed(29)
N <- 10L
B <- 10L
delta <- 0.25

### same right-censored observations
d <- rfree1way(N, delta = c(0, delta), blocks = B, link = "cloglog")
r <- rfree1way(N, delta = c(0, delta), blocks = B, link = "cloglog")
d$y <- Surv(pmin(d$y, r$y), event = d$y < r$y)
## all blocks
ft <- free1way(y ~ groups | blocks, 
               data = d, link = "cloglog")

ex <- expression({

    start <- do.call("c", ft$intercepts)
    fix <- c(0, 0) ### H0
    KM <- ft$profile(c(fix, start), 1:2)$intercepts

    ## compare with Kaplan-Meier
    km <- survfit(y ~ 1, data = subset(d, blocks == paste0("block", b)))
    p1 <- km$surv[km$n.event > 0]
    p0 <- exp(-exp(KM[[b]]))
    ret <- isTRUE(all.equal(p1[seq_along(p0)], p0))

    ## only one block
    ft2 <- free1way(y ~ groups, data = subset(d, blocks == paste0("block", b)),
                    link = "cloglog")
    start <- do.call("c", ft2$intercepts)
    fix <- c(0, 0)
    KM <- ft2$profile(c(fix, start), 1:2)$intercepts
    p2 <- exp(-exp(KM[[1]]))

    ## nonparametric log-likelihoods
    ll1 <- sum(log(diff(1 - c(1, km$surv[km$n.event > 0])))) + 
           sum(log(km$surv[km$n.event == 0]))
    ll2 <- ft2$profile(c(fix, start), 1:2)$value
    ret <- ret && isTRUE(all.equal(ll1 + ll2, 0))

    ret <- ret && isTRUE(all.equal(p1[seq_along(p2)], p2))
})

stopifnot(all(sapply(1:B, function(b) eval(ex))))

