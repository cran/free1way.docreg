
pkgs <- c("free1way.docreg", "brglm2")
if (!all(sapply(pkgs, require, character.only = TRUE, warn.conflicts = FALSE)))
    quit()

### glm fitting and deprecated .structure calls in brglm2
options(digits = 5, warn = -1L)
set.seed(29)

### check log-determinant on simulated data
d <- rfree1way(n = 5, blocks = 3, delta = 1)
ft <- free1way(y ~ groups | blocks, data = d, logdet = TRUE)

H <- get("he", envir = environment(ft$profile))(ft$par)
all.equal(unname(c(ft$logdet)), 
                 c(determinant(as.matrix(H), log = TRUE)$modulus))

### UCBAdmissions (large sample)
ft <- free1way(UCBAdmissions, logdet = TRUE)

H <- get("he", envir = environment(ft$profile))(ft$par)
all.equal(unname(c(ft$logdet)), 
                 c(determinant(as.matrix(H), log = TRUE)$modulus))

ucb <- as.data.frame(UCBAdmissions)

lr <- glm(Admit ~ 0 + Dept + Gender, data = ucb, weights = ucb$Freq, 
          family = binomial(link = "logit"))

log(mantelhaen.test(UCBAdmissions)$estimate)
coef(ft)
rev(coef(lr))[1L]

logLik(ft)
logLik(lr)
diag(vcov(ft))
rev(diag(vcov(lr)))[1L]

confint(ft, test = "Wald")
confint.default(lr)[names(rev(coef(lr))[1L]),]
confint(ft, test = "Rao")
confint(lr, test = "Rao")[names(rev(coef(lr))[1L]),]
confint(ft, test = "LRT")
confint(lr, test = "LRT")[names(rev(coef(lr))[1L]),]

## with penalisation
ftJ <- free1way(UCBAdmissions, MPL_Jeffreys = TRUE, logdet = TRUE)

H <- get("he", envir = environment(ftJ$profile))(ftJ$par)
all.equal(unname(c(ftJ$logdet)), 
                 c(determinant(as.matrix(H), log = TRUE)$modulus))

ucb <- as.data.frame(UCBAdmissions)

Flr <- glm(Admit ~ 0 + Dept + Gender, data = ucb, weights = ucb$Freq, 
           family = binomial(link = "logit"), method = "brglmFit", type = "MPL_Jeffreys")

log(mantelhaen.test(UCBAdmissions)$estimate)
coef(ft)
rev(coef(lr))[1L]

diag(vcov(ft))
rev(diag(vcov(lr)))[1L]

confint(ftJ, test = "Wald")
### brglm2:::confint.brglmFit -> stats:::confint.default
confint(Flr)[names(rev(coef(Flr))[1L]),]

### Rabbits data
example(mantelhaen.test, echo = FALSE)
Rabbits <- aperm(as.table(Rabbits), c(2, 1, 3))

mantelhaen.test(Rabbits, correct = FALSE)

ft <- free1way(Rabbits, logdet = TRUE)

H <- get("he", envir = environment(ft$profile))(ft$par)
all.equal(unname(c(ft$logdet)), 
                 c(determinant(as.matrix(H), log = TRUE)$modulus))

rab <- as.data.frame(Rabbits)

lr <- glm(Response ~ 0 + Penicillin.Level + 
                     Delay, data = rab, weights = rab$Freq, 
          family = binomial(link = "logit"))

coef(ft)
rev(coef(lr))[1L]

logLik(ft)
logLik(lr)

diag(vcov(ft))
rev(diag(vcov(lr)))[1L]

confint(ft, test = "Wald")
confint(ft, test = "Rao")
confint(lr, test = "Rao")[names(rev(coef(lr))[1L]),]
confint(ft, test = "LRT")
confint(lr, test = "LRT")[names(rev(coef(lr))[1L]),]

## with penalisation
ftJ <- free1way(Rabbits, MPL_Jeffreys = TRUE, logdet = TRUE)

H <- get("he", envir = environment(ftJ$profile))(ftJ$par)
all.equal(unname(c(ftJ$logdet)), 
                 c(determinant(as.matrix(H), log = TRUE)$modulus))

Flr <- glm(Response ~ 0 + Penicillin.Level + Delay, data = rab, weights = rab$Freq, 
           family = binomial(link = "logit"), method = "brglmFit", type = "MPL_Jeffreys")

cf <- ftJ$par
cbind(c(-cf[seq_along(cf)[-1L]], cf[1L]), coef(Flr))

log(mantelhaen.test(Rabbits)$estimate)
coef(ftJ)
rev(coef(Flr))[1L]

diag(vcov(ftJ))
rev(diag(vcov(Flr)))[1L]

H <- get("he", envir = environment(ftJ$profile))(ftJ$par)
max(abs(solve(H)[c(seq_along(cf)[-1L], 1), c(seq_along(cf)[-1L], 1)] - vcov(Flr)))

confint(ftJ, test = "Wald")
confint(Flr)[names(rev(coef(Flr))[1L]),]

### repeated measures
d <- rfree1way(n = 1, blocks = 25, delta = 1)
d$y <- factor(d$y < .5)

ft <- free1way(y ~ groups | blocks, data = d, logdet = TRUE)

H <- get("he", envir = environment(ft$profile))(ft$par)
all.equal(unname(c(ft$logdet)), 
                 c(determinant(as.matrix(H), log = TRUE)$modulus))

lr <- glm(y ~ 0 + blocks + groups, data = d, 
          family = binomial(link = "logit"))

coef(ft)
rev(coef(lr))[1L]

logLik(ft)
logLik(lr)

diag(vcov(ft))
rev(diag(vcov(lr)))[1L]

confint(ft, test = "Wald")
confint(ft, test = "Rao")
confint(lr, test = "Rao")[names(rev(coef(lr))[1L]),]
confint(ft, test = "LRT")
confint(lr, test = "LRT")[names(rev(coef(lr))[1L]),]

## with penalisation
ftJ <- free1way(y ~ groups | blocks, data = d, MPL_Jeffreys = TRUE, logdet = TRUE)

H <- get("he", envir = environment(ftJ$profile))(ftJ$par)
all.equal(unname(c(ftJ$logdet)), 
                 c(determinant(as.matrix(H), log = TRUE)$modulus))

Flr <- glm(y ~ 0 + blocks + groups, data = d, 
           family = binomial(link = "logit"), method = "brglmFit", type = "MPL_Jeffreys")

cf <- ftJ$par
cbind(c(-cf[seq_along(cf)[-1L]], cf[1L]), coef(Flr))

log(mantelhaen.test(xtabs(~ y + groups + blocks, data = d))$estimate)
coef(ftJ)
rev(coef(Flr))[1L]

diag(vcov(ftJ))
rev(diag(vcov(Flr)))[1L]

H <- get("he", envir = environment(ftJ$profile))(ftJ$par)
max(abs(solve(H)[c(seq_along(cf)[-1L], 1), c(seq_along(cf)[-1L], 1)] - vcov(Flr)))

confint(ftJ, test = "Wald")
confint(Flr)[names(rev(coef(Flr))[1L]),]

### with powers of Jeffrey's prior
ftJ <- free1way(y ~ groups | blocks, data = d, MPL_Jeffreys = 2 * 1.5, logdet = TRUE)

H <- get("he", envir = environment(ftJ$profile))(ftJ$par)
all.equal(unname(c(ftJ$logdet)), 
                 c(determinant(as.matrix(H), log = TRUE)$modulus))

Flr <- glm(y ~ 0 + blocks + groups, data = d, 
           family = binomial(link = "logit"), method = "brglmFit", 
           type = "MPL_Jeffreys", a = 1.5)

cf <- ftJ$par
cbind(c(-cf[seq_along(cf)[-1L]], cf[1L]), coef(Flr))

log(mantelhaen.test(xtabs(~ y + groups + blocks, data = d))$estimate)
coef(ftJ)
rev(coef(Flr))[1L]

diag(vcov(ftJ))
rev(diag(vcov(Flr)))[1L]

H <- get("he", envir = environment(ftJ$profile))(ftJ$par)
max(abs(solve(H)[c(seq_along(cf)[-1L], 1), c(seq_along(cf)[-1L], 1)] - vcov(Flr)))

confint(ftJ, test = "Wald")
### brglm2:::confint.brglmFit -> stats:::confint.default
confint(Flr)[names(rev(coef(Flr))[1L]),]
