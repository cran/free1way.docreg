
# scale tests

scaleTest <- function(object, B = 0, ...) {

    alternative <- "two.sided"
    y <- object$table
    d <- dim(y)
    ret <- object
    cf <- ret$par
    cf[idx <- seq_len(d[2L] - 1L)] <- 0
    pr <- ret$profile(cf, idx)
    res <- - pr$negscaleresiduals

    # Strasser Weber
    
    .SW <- function(res, xt) 
    {

        if (length(dim(xt)) == 3L) {
            res <- matrix(res, nrow = dim(xt)[1L], ncol = dim(xt)[3])
            STAT <-  Exp <- Cov <- 0
            for (b in seq_len(dim(xt)[3L])) {
                sw <- .SW(res[,b, drop = TRUE], xt[,,b, drop = TRUE])
                STAT <- STAT + sw$Statistic
                Exp <- Exp + sw$Expectation
                Cov <- Cov + sw$Covariance
            }
            return(list(Statistic = STAT, Expectation = as.vector(Exp),
                        Covariance = Cov))
        }

        Y <- matrix(res, ncol = 1, nrow = length(xt))
        weights <- c(xt)
        x <- gl(ncol(xt), nrow(xt))
        X <- model.matrix(~ x, data = data.frame(x = x))[,-1L,drop = FALSE]

        w. <- sum(weights)
        wX <- weights * X
        wY <- weights * Y
        ExpX <- colSums(wX)
        ExpY <- colSums(wY) / w.
        CovX <- crossprod(X, wX)
        Yc <- t(t(Y) - ExpY)
        CovY <- crossprod(Yc, weights * Yc) / w.
        Exp <- kronecker(ExpY, ExpX)
        Cov <- w. / (w. - 1) * kronecker(CovY, CovX) -
               1 / (w. - 1) * kronecker(CovY, tcrossprod(ExpX))
        STAT <- crossprod(X, wY)
        list(Statistic = STAT, Expectation = as.vector(Exp),
             Covariance = Cov)
    }
    
    # resampling
    
    .resample <- function(res, xt, B = 10000) 
    {

        if (length(dim(xt)) == 2L)
            xt <- as.table(array(xt, dim = c(dim(xt), 1)))

        res <- matrix(res, nrow = dim(xt)[1L], ncol = dim(xt)[3L])
        stat <- 0
        ret <- .SW(res, xt)
        if (dim(xt)[2L] == 2L) {
            ret$testStat <- c((ret$Statistic - ret$Expectation) / 
                              sqrt(c(ret$Covariance)))
        } else {
            ES <- ret$Statistic - ret$Expectation
            ret$testStat <- sum(ES * solve(ret$Covariance, ES))
        }
        ret$DF <- dim(xt)[2L] - 1L

        if (B) {
            for (j in 1:dim(xt)[3L]) {
               rt <- r2dtable(B, r = rowSums(xt[,,j]), c = colSums(xt[,,j]))
               stat <- stat + vapply(rt, 
                   function(x) .colSums(x[,-1L, drop = FALSE] * res[,j], 
                                        m = nrow(x), n = ncol(x) - 1L), 
                                     FUN.VALUE = rep(0, dim(xt)[[2L]] - 1L))
            }
            if (dim(xt)[2L] == 2L) {
                 ret$permStat <- (stat - ret$Expectation) / 
                                  sqrt(c(ret$Covariance))
            } else {
                ES <- matrix(stat, ncol = B) - ret$Expectation
                ret$permStat <- .colSums(ES * solve(ret$Covariance, ES), 
                                         m = dim(xt)[[2L]] - 1L, n = B)
            }
        }
        ret
    }
    

    if (length(dim(y)) == 3L) y <- y[,,ret$strata, drop = FALSE]
    if (length(dim(y)) == 4L) {
        y <- y[,,ret$strata,, drop = FALSE]
        dy <- dim(y)
        dy[1] <- dy[1] * 2
        y <- apply(y, 3, function(x) rbind(x[,,"TRUE"], x[,,"FALSE"]), simplify = FALSE)
        y <- array(unlist(y), dim = dy[1:3])
    }
    ret$perm <- .resample(res, y, B = B)

    STATISTIC <- ret$perm$testStat
    x <- ret
    cf <- coef(object)
    sc <- 0
    # Permutation p-values
    
    if (!is.null(x$exact)) {
        PVAL <- switch(alternative,
                       "two.sided" = pmin(1, 2 * min(c(x$exact$ple(sc), 
                                                       x$exact$pgr(sc)))),
                       "less" = x$exact$ple(sc),
                       "greater" = x$exact$pgr(sc))
    } else {
        .pm <- function(x) sum(x) / length(x) 
        ps <- x$perm$permStat

        .GE <- function(x, y)
            (y - x) <= sqrt(.Machine$double.eps)

        .LE <- function(x, y)
            (x - y) <= sqrt(.Machine$double.eps)

        if (alternative == "two.sided" && length(cf) > 1L) {
            if (!is.null(ps)) {
                PVAL <- .pm(.GE(ps, STATISTIC))
            } else {
                DF <- c("df" = x$perm$DF)
                PVAL <- pchisq(STATISTIC, df = DF, lower.tail = FALSE)
            }
        } else {
            if (!is.null(ps)) {
                PVALle <- .pm(.LE(ps, STATISTIC))
                PVALge <- .pm(.GE(ps, STATISTIC))
                if (alternative == "two.sided")
                    PVAL <- 2 * min(c(PVALle, PVALge))
                else if (alternative == "less")
                    PVAL <- PVALle
                else
                    PVAL <- PVALge
            } else {
                if (alternative == "two.sided")
                    PVAL <- pchisq(STATISTIC^2, df = 1, lower.tail = FALSE)
                else
                    PVAL <- pnorm(STATISTIC, 
                                  lower.tail = alternative == "less")
            }
        }
    }
    

    RVAL <- list(statistic = STATISTIC, parameter = NULL, p.value = PVAL, 
        null.value = 1, alternative = alternative, method = "Permutation scale test",
        data.name = x$data.name)
    class(RVAL) <- "htest"
    return(RVAL)
}

