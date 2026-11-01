
library("free1way.docreg", warn.conflicts = FALSE)
set.seed(290875)

### bugs fixed in 4.6.1, so test with R 4.6.1 or later only
if (compareVersion(paste0(version$major, ".", version$minor), "4.6.1") >= 0)
{
### strata_ratio didn't work properly in free1way.docreg 1.0-0 / R 4.6.0
### fixed 2026-06-05
df <- rfree1way(10, delta = .5, strata_ratio = c(2, 3), blocks = 3)
(xt <- xtabs(~ groups + blocks, data = df))
stopifnot(isTRUE(max(abs(xt - 10 * matrix(rep(1:3, 2), byrow = TRUE, nrow = 2))) == 0))

### while we are at it, test alloc_ratio together with strata_ratio
df <- rfree1way(10, delta = .5, strata_ratio = c(2, 3), alloc_ratio = 2, blocks = 3)
(xt <- xtabs(~ groups + blocks, data = df))
stopifnot(isTRUE(max(abs(xt - 10 * matrix(rep(1:3, 2), byrow = TRUE, nrow = 2) * c(1, 2))) == 0))

### prob arg with zero prob categories resulted in 
### Error in cut.default(runif(colsums[k]), breaks = c(-Inf, p[, k])) : 
#  'breaks' are not unique
power.free1way.test(n = 10, prob = matrix(c(0, 0, 1, 2, 3), ncol = 1),
                    delta = 1)
}

### bugs fixed in 4.6.2, so test with R 4.6.2 or later only
if (compareVersion(paste0(version$major, ".", version$minor), "4.6.2") >= 0)
{
### Overlap coefficients are monotone increasing for negative arguments and
### decreasing for positive ones.  Confidence intervals for shift terms
### including zero should include 1 for the overlap coefficient
d <- rfree1way(10, delta = c(0, -2, 2))
fw <- free1way(y ~ groups, data = d)
CI <- confint(fw, test = "Wald", what = "OVL")
stopifnot(all(CI[,1] < CI[,2]))
}
