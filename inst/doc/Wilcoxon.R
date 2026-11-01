### R code from vignette source 'Wilcoxon.Rnw'

###################################################
### code chunk number 1: packages
###################################################
pkgs <- c("free1way.docreg", "ordinalCont", "lattice", "latticeExtra", 
          "MASS", "memisc", "survival", "bibtex")
OK <- sapply(pkgs, require, character.only = TRUE)
if (!all(OK)) {
    install.packages(pkgs[!OK])
    OK <- sapply(pkgs, require, character.only = TRUE)
    if (!all(OK))
        stop("package(s) ", paste(pkgs[!OK], collapse = ", "), " not available")
}
set.seed(290875)
options(prompt = "R> ", continue = "+  ", width = 70, 
        useFancyQuotes = FALSE, digits = 6)


###################################################
### code chunk number 2: tools
###################################################
## some helper functions
frmt1 <- function(x, digits = 1) {
    if (!is.numeric(x)) return(x)
    formatC(round(x, digits = digits), digits = digits, format = "f") 
}
frmt3 <- function(x) 
    frmt1(x, digits = 3)
frmtci <- function(x, digits = 3) {
    if (!is.numeric(x)) return(x)
    if (length(x) != 2) stop("not a confidence interval")
    return(paste("(", frmt1(x[1], digits = digits), 
                 ",", frmt1(x[2], digits = digits), ")"))
}


###################################################
### code chunk number 3: pain
###################################################
pain_df <- neck_pain
idt <- xtabs(~ id, data = pain_df)
miss <- names(idt)[idt < 3]
#pain_df <- subset(pain_df, !id %in% miss)
#pain_df$id <- factor(pain_df$id)
levels(pain_df$laser) <- c("Placebo", "Active")
levels(pain_df$time) <- c("Baseline", "7 weeks", "12 weeks")
pain <- rbind(subset(pain_df, laser == levels(pain_df$laser)[1]),
              subset(pain_df, laser == levels(pain_df$laser)[2]))
ties12 <- sum(duplicated(subset(pain, time == "12 weeks")$vas))
pain$vasB <- cut(pain$vas, breaks = c(-Inf, .6, Inf), 
                 ordered_result = FALSE, 
                 labels = c("None to moderate", "(Moderately) severe"))
cols <- trellis.par.get("superpose.line")$col[1:2]
names(cols) <- levels(pain$laser)


###################################################
### code chunk number 4: Bauer
###################################################
Bauer <- matrix(c(
2025,	63,
2024,	58,
2023,	54,
2022,	60,
2021,	50,
2020,	34,
2019,	41,
2018,	41,
2017,	32,
2016,	33,
2015,	20,
2014,	19,
2013,	16,
2012,	11,
2011,	16,
2010,	11,
2009,	14,
2008,	2,
2007,	2,
2006,	2,
2005,	2,
2003,	1,
2002,	1,
2001,	1,
1997,	3,
1996,	2,
1995,	1,
1992,	2,
1991,	2,
1990,	1,
1989,	1,
1988,	2,
1985,	1,
1982,	1,
1981,	1,
1980,	1,
1976,	3,
1974,	1,
1973,	1), byrow = TRUE, ncol = 2)
Bauer <- Bauer[order(Bauer[,1]),]
yrs <- Bauer[,1]
x <- Bauer[,2]
pre2000 <- sum(x[yrs < 2000])
plot(yrs, cumsum(x), type = "s", xlab = "Year Published",
     ylab = "Cumulative Number of Citations", las = 1, xlim = c(1970, max(yrs) + 1))
arrows(2000.5, sum(x) / 3, 2000.5, 28, length = .15)
text(2000.5, sum(x) / 3 + 55, "R 1.1.0 implements 
wilcox.test(,conf.int = TRUE)")
arrows(1972, sum(x) / 3, 1972, 5, length = .15)
text(1972, sum(x) / 3, "Publication", pos = 4)


###################################################
### code chunk number 5: samplesize
###################################################
toLatex(ftable(xtabs(~ laser + time, data = pain)))


###################################################
### code chunk number 6: plot
###################################################
print(xyplot(vas ~ time | laser, data = pain, 
       groups = id, type = "b", col = rgb(.1, .1, .1, .1),
       lwd = 2, layout = c(2, 1),
       ylab = "Neck pain (on visual analogue scale)", 
       xlab = "Examinations"))


###################################################
### code chunk number 7: ecdf
###################################################
print(ecdfplot(~ vas | time, data = pain, groups = laser,
         xlab = "Neck pain (on visual analogue scale)", 
         panel = function(...) {
             panel.ecdfplot(...)
             panel.abline(a = 0, b = 1, col = "grey")
         },
         auto.key = list(columns = 2, position = "top")))


###################################################
### code chunk number 8: wilcox.test
###################################################
(w <- wilcox.test(vas ~ laser, data = pain, subset = time == "12 weeks",
                  exact = TRUE, conf.int = TRUE))


###################################################
### code chunk number 9: polr
###################################################
w12 <- subset(pain, time == "12 weeks")
w12$vasE <- ordered(w12$vas)
p <- polr(vasE ~ laser, data = w12, Hess = TRUE)
summary(p)$coefficients["laserActive",]
confint(p)


###################################################
### code chunk number 10: wpolr
###################################################
wilcox.test(vas ~ laser, data = pain, subset = time == "12 weeks",
            exact = FALSE)$p.value


###################################################
### code chunk number 11: wpolrci
###################################################
confint(p, level = 1 - w$p.value)


###################################################
### code chunk number 12: po
###################################################
w12 <- cbind(subset(pain, time == "12 weeks"), ev = TRUE)
plot(survfit(Surv(vas, ev) ~ laser, data = w12),
     conf.int = FALSE, mark.time = FALSE, 
     fun = function(p) qlogis(p, lower.tail = FALSE),
     col = cols, las = 1, main = "12 weeks",
     xlab = "Neck pain (on visual analogue scale)",
     ylab = "Distribution function on logit scale")
legend("bottomright", legend = levels(w12$laser), col = cols,
       lty = 1, bty = "n", title = "Empirical")


###################################################
### code chunk number 13: free1waytmp
###################################################
f <- free1way(vas ~ laser, data = pain, subset = time == "12 weeks")


###################################################
### code chunk number 14: free1way
###################################################
f <- free1way(vas ~ laser, data = pain, subset = time == "12 weeks", 
              exact = TRUE)


###################################################
### code chunk number 15: loglik
###################################################
### polr
logLik(p)
### free1way
logLik(f)
### number of log-odds ratio + intercept parameters
length(coef(f)) + length(f$intercepts[[1L]])


###################################################
### code chunk number 16: coef
###################################################
all.equal(coef(p), coef(f))
all.equal(p$zeta, f$intercepts[[1L]], check.attributes = FALSE)
all.equal(vcov(p)[names(coef(p)),names(coef(p)), drop = FALSE], 
          vcov(f), check.attributes = FALSE)


###################################################
### code chunk number 17: fplot
###################################################
w12 <- cbind(subset(pain, time == "12 weeks"), ev = TRUE)
plot(survfit(Surv(vas, ev) ~ laser, data = w12),
     conf.int = FALSE, mark.time = FALSE, 
     fun = function(p) qlogis(p, lower.tail = FALSE),
     col = cols, las = 1, main = "12 weeks",
     xlab = "Neck pain (on visual analogue scale)",
     ylab = "Distribution function on logit scale")
legend("bottomright", legend = levels(w12$laser), col = cols,
       lty = 1, bty = "n", title = "Empirical")
su <- sort(unique(w12$vas))
lines(su[-length(su)], f$intercepts[[1]], 
      type = "s", lty = 2, col = cols["Placebo"])
lines(su[-length(su)], f$intercepts[[1]] - coef(f), 
      type = "s", lty = 2, col = cols["Active"])
legend("topleft", legend = levels(w12$laser), col = cols, 
       lty = 2, bty = "n", title = "Model-based")
suq <- quantile(w12$vas, prob = 1:3 / 4)
out <- sapply(suq, function(q) {
    px <- f$intercepts[[1]][which.min((su - q)^2)]
    arrows(q, px, q, px - coef(f), code = 3, length = .1, angle = 45)
    text(q, px - coef(f) + .5, labels = expression(hat(delta)))
})


###################################################
### code chunk number 18: free1wayperm
###################################################
(sf <- summary(f, test = "Permutation"))
coef(f)
confint(f, test = "Permutation")


###################################################
### code chunk number 19: free1waypermOR
###################################################
(OR <- exp(coef(f)))
(ciOR <- exp(confint(f, test = "Permutation")))


###################################################
### code chunk number 20: free1waypermAUC
###################################################
coef(f, what = "AUC")
confint(f, test = "Permutation", what = "AUC")


###################################################
### code chunk number 21: gof
###################################################
x <- subset(pain, subset = time == "12 weeks")
with(split(x$vas, x$laser), 
     ppplot(Active, Placebo, conf.level = .95, 
            main = "Visual Analogue Scale at 12 weeks"))


###################################################
### code chunk number 22: free1waypermOVL
###################################################
coef(f, what = "OVL")
confint(f, test = "Permutation", what = "OVL")


###################################################
### code chunk number 23: free1waypermCI
###################################################
confint(f, level = 1 - sf$p.value, test = "Permutation")


###################################################
### code chunk number 24: free1wayRao
###################################################
summary(f, test = "Rao")
confint(f, test = "Rao")


###################################################
### code chunk number 25: free1wayLRT
###################################################
summary(f, test = "LRT")
confint(f, test = "LRT")


###################################################
### code chunk number 26: free1wayWald
###################################################
summary(f, test = "Wald")
confint(f, test = "Wald")


###################################################
### code chunk number 27: power
###################################################
(pt <- power.free1way.test(delta = coef(f), power = .8, sig.level = .05))


###################################################
### code chunk number 28: simulation (eval = FALSE)
###################################################
## ex <- expression({
##     df <- rfree1way(n = pt$n, delta = pt[["log-odds ratio"]])
##     summary(free1way(y ~ groups, data = df), 
##             test = "Permutation")$p.value
## })
## pvals <- replicate(1000, eval(ex))
## binom.test(x = sum(pvals < pt[["sig.level"]]), n = length(pvals), 
##            p = pt[["power"]])$conf.int


###################################################
### code chunk number 29: vasO
###################################################
pain$vasO <- cut(pain$vas, breaks = c(-Inf, .3, .6, .8, Inf), 
                 ordered_result = TRUE, 
                 labels = c("None to mild", "Moderate", 
                            "Moderately severe", "Severe"))


###################################################
### code chunk number 30: vasBL
###################################################
xB <- subset(pain, time == "Baseline")
pain$vasBL <- xB$vasO[match(pain$id, xB$id)]
class(pain$vasBL) <- "factor"


###################################################
### code chunk number 31: free1wayS
###################################################
fb <- free1way(vas ~ laser | vasBL, data = pain, 
               subset = time == "12 weeks", B = 10000)
summary(fb, test = "Permutation")
exp(coef(fb))
exp(confint(fb, test = "Permutation"))


###################################################
### code chunk number 32: free1wayO
###################################################
fO <- free1way(vasO ~ laser | vasBL, data = pain, 
               subset = time == "12 weeks", B = 10000)
summary(fO, test = "Permutation")
exp(coef(fO))
exp(confint(fO, test = "Permutation"))


###################################################
### code chunk number 33: placebo
###################################################
(btab <- xtabs(~ vasO + vasBL, data = pain, 
               subset = laser == "Placebo" & time == "12 weeks")[,4:1])
colSums(btab)
(ptab <- prop.table(btab, margin = 2))


###################################################
### code chunk number 34: powerb
###################################################
(pbt <- power.free1way.test(prob = ptab[,-1], strata_ratio = c(2, 2),
                            delta = coef(fO), power = .8, 
                            sig.level = .05))


###################################################
### code chunk number 35: mh
###################################################
mh <- free1way(vasB ~ laser | vasBL, data = pain, 
               subset = time == "12 weeks")
summary(mh, test = "Permutation")$p.value
xt <- xtabs(~ vasB + laser + vasBL, 
            data = pain, subset = time == "12 weeks")
mantelhaen.test(xt, correct = FALSE, exact = FALSE)$p.value


###################################################
### code chunk number 36: mhOR
###################################################
### Mantel-Haenszel
mantelhaen.test(xt, correct = FALSE, exact = FALSE)$estimate
### Conditional ML
mantelhaen.test(xt, correct = FALSE, exact = TRUE)$estimate
### ML
exp(coef(mh))


###################################################
### code chunk number 37: indivplot
###################################################
pain$vasBLorg <- xB$vas[match(pain$id, xB$id)]
print(xyplot(vas ~ vasBLorg | time + laser, 
       data = subset(pain, time != "Baseline"),
       panel = function(x, y, ...) {
           panel.abline(a = 0, b = 1, col = "grey")
           panel.xyplot(x, y, ...)
           ltext(.1, .9, label = paste0(round(mean(y < x) * 100), "%"))
       }, xlab = "Baseline neck pain (on visual analogue scale)",
       ylab = "Neck pain (on visual analogue scale)"))


###################################################
### code chunk number 38: friedmanP
###################################################
fP <- free1way(vas ~ time | id, data = pain, 
               subset = laser == "Placebo" & !(id %in% miss))
unlist(summary(fP, test = "Permutation")[c("statistic", "p.value")])
painP <- subset(pain, laser == "Placebo" & !(id %in% miss))
painP$id <- painP$id[,drop = TRUE]
unlist(friedman.test(vas ~ time | id, 
                     data = painP)[c("statistic", "p.value")])


###################################################
### code chunk number 39: mcnemar
###################################################
summary(free1way(vasB ~ time | id, data = pain, 
                 subset = laser == "Placebo" & time != "7 weeks"), 
        test = "Permutation")$p.value
xt <- xtabs(~ vasBL + vasB, data = pain, 
             subset = laser == "Placebo" & time != "7 weeks")
xt <- rbind(colSums(xt[1:2,]), colSums(xt[3:4,]))
mcnemar.test(xt, correct = FALSE)$p.value


###################################################
### code chunk number 40: friedmanA
###################################################
fA <- free1way(vas ~ time | id, data = pain, 
               subset = laser == "Active" & !(id %in% miss))
unlist(summary(fA, test = "Permutation")[c("statistic", "p.value")])
painA <- subset(pain, laser == "Active" & !(id %in% miss))
painA$id <- painA$id[,drop = TRUE]
unlist(friedman.test(vas ~ time | id, 
                     data = painA)[c("statistic", "p.value")])


###################################################
### code chunk number 41: MPL-A (eval = FALSE)
###################################################
## fA <- free1way(vas ~ time | id, data = pain, 
##                subset = laser == "Active" & !(id %in% miss))
## fAp <- free1way(vas ~ time | id, data = pain, 
##                 subset = laser == "Active" & !(id %in% miss), 
##                 MPL_Jeffreys = TRUE)
## cbind(unpenalised = coef(fA), penalised = coef(fAp))


###################################################
### code chunk number 42: MPL-P (eval = FALSE)
###################################################
## fP <- free1way(vas ~ time | id, data = pain, 
##                subset = laser == "Placebo" & !(id %in% miss))
## fPp <- free1way(vas ~ time | id, data = pain, 
##                 subset = laser == "Placebo" & !(id %in% miss), 
##                 MPL_Jeffreys = TRUE)
## cbind(unpenalised = coef(fP), penalised = coef(fPp))


###################################################
### code chunk number 43: MPL-diff (eval = FALSE)
###################################################
## cbind(unpenalised = coef(fA) - coef(fP), 
##       penalised = coef(fAp) - coef(fPp))


###################################################
### code chunk number 44: bib
###################################################
if (file.exists("REFERENCES.bib"))
    file.remove("REFERENCES.bib")
pkgs <- c("free1way.docreg", "tram")
bib <- do.call("c", lapply(pkgs, function(pkg) read.bib(package = pkg)))
bib <- bib[!duplicated(names(bib))]
if ("R" %in% names(bib))
    bib[["R"]] <- NULL
write.bib(c("rms", "tram", "MASS", "free1way.docreg", "Matrix", "survival"))
bib["Harrell2015RMS"]$doi <- gsub("_", "\\\\_", bib["Harrell2015RMS"]$doi)
write.bib(bib, file = "REFERENCES.bib")
R <- unclass(toBibtex(citation()))
R[[1]] <- "@Manual{R,"
cat(R, file = "Rpackages.bib", append = TRUE, sep = "\n")


###################################################
### code chunk number 45: vign (eval = FALSE)
###################################################
## vignette("free1way", package = "free1way.docreg")


###################################################
### code chunk number 46: Jeffrey-nlminb (eval = FALSE)
###################################################
## fPp <- free1way(vas ~ time | id, data = pain, 
##                 subset = laser == "Placebo", 
##                 MPL_Jeffreys = TRUE, logdet = TRUE)


###################################################
### code chunk number 47: Jeffrey-chk (eval = FALSE)
###################################################
## pll <- function(theta) {
##     ret <- fPp$profile(theta, fix = seq_along(theta))
##     ret$value - .5 * c(ret$logdet)
## }
## logLik(fPp)
## -pll(fPp$par)


###################################################
### code chunk number 48: Jeffrey-optim (eval = FALSE)
###################################################
## op <- optim(fPp$par, fn = pll)
## logLik(fPp)
## -op$value
## coef(fPp)
## op$par[1:2]


###################################################
### code chunk number 49: Jeffrey-dooptim (eval = FALSE)
###################################################
## optNM <- function(start, objective, control = list(trace = FALSE), ...) {
##     ctrl <- list(trace = control$trace)
##     ret <- optim(start, fn = objective, control = ctrl)
##     ret$objective <- ret$value
##     ret
## }
## fPpNM <- free1way(vas ~ time | id, data = pain, dooptim = "optNM",
##     subset = laser == "Placebo", 
##     control = list(optNM = list(trace = FALSE, iter.max = 100)),
##     MPL_Jeffreys = TRUE, logdet = TRUE)
## logLik(fPp)
## logLik(fPpNM)


###################################################
### code chunk number 50: binJeffrey (eval = FALSE)
###################################################
## library("brglm2")
## xt <- xtabs(~ vasB + time + id, data = pain, 
##             subset = laser == "Placebo" & time != "7 weeks",
##             drop.unused.levels = TRUE)
## dxt <- as.data.frame(xt)
## Firth <- glm(vasB ~ 0 + id + time, 
##              data = dxt, weights = dxt$Freq, 
##              family = binomial(link = "logit"), 
##              method = "brglmFit", type = "MPL_Jeffreys")
## exp(coef(Firth)["time12 weeks"])
## exp(coef(fJ <- free1way(xt, MPL_Jeffreys = TRUE)))
## mantelhaen.test(xt)$estimate
## ### biased ML, the default
## exp(coef(free1way(xt, MPL_Jeffreys = FALSE)))


###################################################
### code chunk number 51: Firth (eval = FALSE)
###################################################
## max(abs(coef(Firth)[-length(coef(Firth))] + fJ$par[-1L]))


###################################################
### code chunk number 52: small-wt
###################################################
N <- 5
d <- data.frame(groups = gl(2, N))
d$y <- seq_len(2 * N)
wilcox.test(y ~ groups, data = d, exact = TRUE)$p.value


###################################################
### code chunk number 53: small-f1w
###################################################
f <- free1way(y ~ groups, data = d, exact = TRUE)
summary(f, test = "Permutation")$p.value


###################################################
### code chunk number 54: small-f1w-he
###################################################
f$hessian


###################################################
### code chunk number 55: small-f1w-ll
###################################################
coef(f)
theta <- f$par[-1L]
delta <- 20:60
logLik(f)
range(sapply(delta, function(d)
    -f$profile(c(d, theta), fix = 1L)$value))


###################################################
### code chunk number 56: addgroups
###################################################
d <- rfree1way(25, delta = c(.5, .75))

### identical
cbind(coef(lm(y ~ groups, data = d))[-1L],
      c(coef(lm(y ~ groups, data = d, subset = groups != "C"))[-1L],
        coef(lm(y ~ groups, data = d, subset = groups != "B"))[-1L]))

### not identical
cbind(coef(free1way(y ~ groups, data = d)),
      c(coef(free1way(y ~ groups, data = d, subset = groups != "C")),
        coef(free1way(y ~ groups, data = d, subset = groups != "B"))))


###################################################
### code chunk number 57: sessionInfo
###################################################
toLatex(sessionInfo())
