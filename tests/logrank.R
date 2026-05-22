
set.seed(29)
options(digits = 5)

pkgs <- c("free1way.docreg", "survival", "tram")
if (!all(sapply(pkgs, require, character.only = TRUE)))
    quit()

beta <- 1
N <- 10
w <- gl(2, N)
d <- data.frame(time = rexp(length(w), rate = c(1, 1 + beta)[w]), w = w)
d$cens <- rbinom(length(w), size = 1L, prob = .3) == 1

m <- Polr(R(Surv(time, cens), as.R.ordered = TRUE) ~ w, data = d, method = "cloglog")
coef(m)

vcov(m)
residuals(m)

(cf <- coef(as.mlt(m)))

ft <- with(d, free1way(y = time, groups = w, 
                       event = cens, 
                       link = "cloglog"))
logLik(ft)
logLik(m)

coef(ft)
coef(m)

vcov(ft)
vcov(m)

pr <- ft$par
cumsum(pr[-1])
cf

ft
perm_test(m)

data("GBSG2", package = "TH.data")

survdiff(Surv(time, cens) ~ horTh, data = GBSG2)

m <- coxph(Surv(time, cens) ~ horTh, data = GBSG2)
coef(m)
vcov(m)

perm_test(m)

ft <- with(GBSG2, free1way(Surv(time, cens) ~ horTh,
                           link = "cloglog"))
coef(ft)
vcov(ft)

ft

survdiff(Surv(time, cens) ~ horTh + strata(tgrade), data = GBSG2)

m <- coxph(Surv(time, cens) ~ horTh + strata(tgrade), data = GBSG2)
coef(m)
vcov(m)

# perm_test(m)

ft <- with(GBSG2, free1way(Surv(time, cens) ~ horTh | tgrade, 
                           link = "cloglog"))
coef(ft)
vcov(ft)

ft

