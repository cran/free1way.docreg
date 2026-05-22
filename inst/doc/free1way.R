### R code from vignette source 'free1way.Rnw'

###################################################
### code chunk number 1: citation
###################################################
yr <- format(dt <- as.Date(packageDescription("free1way.docreg")$Date), "%Y")
vs <- packageDescription("free1way.docreg")$Version
title <- "Semiparametrically Efficient Population and Permutation Inference in 
       Distribution-free Stratified $K$-sample Oneway Layouts"
DOI <- paste0("10.32614/CRAN.package.", packageDescription("free1way.docreg")$Package)


###################################################
### code chunk number 2: localfun
###################################################
Nsim <- 100
options(digits = 5)
.table2list <- function(x) 
{

    
dx <- dim(x)
if (length(dx) == 1L)
    stop("incorrect dimensions")
if (length(dx) == 2L)
    x <- as.table(array(x, dim = c(dx, 1)))
dx <- dim(x)
if (length(dx) < 3L)
    stop("incorrect dimensions")
C <- dim(x)[1L]
K <- dim(x)[2L]
B <- dim(x)[3L]
if (C < 2L)
    stop("at least two response categories required")
if (K < 2L)
    stop("at least two groups required")
xrc <- NULL
if (length(dx) == 4L) {
    if (dx[4] == 2L) {
        xrc <- array(x[,,,"FALSE", drop = TRUE], dim = dx[1:3])
        x <- array(x[,,,"TRUE", drop = TRUE], dim = dx[1:3])
    } else {
        stop(gettextf("%s currently only allows independent right-censoring",
                      "free1way"),
             domain = NA)
    }
}


xlist <- xrclist <- vector(mode = "list", length = B)

for (b in seq_len(B)) {
    xb <- matrix(x[,,b, drop = TRUE], ncol = K)
    xw <- rowSums(abs(xb)) > 0
    if (sum(xw) > 1L) {
        ### do not remove last parameter if there are corresponding
        ### right-censored observations
        wm <- which(xw)[sum(xw)]
        if (!is.null(xrc) && any(xrc[wm:dx[1],,b,drop = TRUE] > 0))
            xw[length(xw)] <- TRUE
        xlist[[b]] <- xb[xw,,drop = FALSE]
        Cidx <- rep.int(1L, times = C)
        Cidx[xw] <- Cidx[xw] + seq_len(sum(xw))
        attr(xlist[[b]], "idx") <- Cidx
        if (!is.null(xrc)) {
            ### count right-censored observations between distinct event
            ### times
            cs <- apply(xrc[,,b,drop = TRUE] * (!xw), 2, function(x) 
                diff(c(0, cumsum(x)[xw])))
            xrclist[[b]] <- matrix(xrc[xw,,b,drop = TRUE], ncol = K) + cs
            idx <- seq_len(C)[xw]
            idx <- rep(seq_len(sum(xw)), times = c(idx[1], diff(idx)))
            Cidx <- rep.int(1L, times = C)
            Cidx[seq_along(idx)] <- Cidx[seq_along(idx)] + idx
            attr(xrclist[[b]], "idx") <- Cidx
        }
    }
}
### remove empty blocks
strata <- !vapply(xlist, is.null, NA)
xlist <- xlist[strata]
xrclist <- xrclist[strata]



    ret <- list(xlist = xlist)
    if (!is.null(xrc))
        ret$xrclist <- xrclist
    ret$strata <- strata
    ret
}

.nll <- function(parm, x, mu = 0, rightcensored = FALSE) 
{

    
    bidx <- seq_len(ncol(x) - 1L)
    delta <- c(0, mu + parm[bidx])
    intercepts <- c(-Inf, parm[- bidx], Inf)
    tmb <- intercepts - matrix(delta, nrow = length(intercepts),  
                                      ncol = ncol(x),
                                      byrow = TRUE)
    Ftmb <- F(tmb)
    if (rightcensored) {
        prb <- 1 - Ftmb[- nrow(Ftmb), , drop = FALSE]
    } else {
        prb <- Ftmb[- 1L, , drop = FALSE] - 
               Ftmb[- nrow(Ftmb), , drop = FALSE]
    } 
    
    if (any(prb < .Machine$double.eps^10)) 
        return(Inf)
    return(- sum(x * log(prb)))
}


.nsc <- function(parm, x, mu = 0, rightcensored = FALSE) 
{

    
    bidx <- seq_len(ncol(x) - 1L)
    delta <- c(0, mu + parm[bidx])
    intercepts <- c(-Inf, parm[- bidx], Inf)
    tmb <- intercepts - matrix(delta, nrow = length(intercepts),  
                                      ncol = ncol(x),
                                      byrow = TRUE)
    Ftmb <- F(tmb)
    if (rightcensored) {
        prb <- 1 - Ftmb[- nrow(Ftmb), , drop = FALSE]
    } else {
        prb <- Ftmb[- 1L, , drop = FALSE] - 
               Ftmb[- nrow(Ftmb), , drop = FALSE]
    } 
    

    
    ftmb <- f(tmb)
    zu <- x * ftmb[- 1, , drop = FALSE] / prb
    if (rightcensored) zu[] <- 0 ### derivative of a constant
    zl <- x * ftmb[- nrow(ftmb), , drop = FALSE] / prb
    

    ret <- numeric(length(parm))
    ret[bidx] <- .colSums(zl, m = nrow(zl), n = ncol(zl))[-1L] -
                 .colSums(zu[-nrow(zu),,drop = FALSE], 
                          m = nrow(zu) - 1L, n = ncol(zu))[-1L]
    ret[- bidx] <- .rowSums(zu[-nrow(zu),,drop = FALSE] - 
                            zl[-1,,drop = FALSE], 
                            m = nrow(zu) - 1L, n = ncol(zu))
    return(- ret)
}


.nsr <- function(parm, x, mu = 0, rightcensored = FALSE) 
{

    
    bidx <- seq_len(ncol(x) - 1L)
    delta <- c(0, mu + parm[bidx])
    intercepts <- c(-Inf, parm[- bidx], Inf)
    tmb <- intercepts - matrix(delta, nrow = length(intercepts),  
                                      ncol = ncol(x),
                                      byrow = TRUE)
    Ftmb <- F(tmb)
    if (rightcensored) {
        prb <- 1 - Ftmb[- nrow(Ftmb), , drop = FALSE]
    } else {
        prb <- Ftmb[- 1L, , drop = FALSE] - 
               Ftmb[- nrow(Ftmb), , drop = FALSE]
    } 
    

    
    ftmb <- f(tmb)
    zu <- x * ftmb[- 1, , drop = FALSE] / prb
    if (rightcensored) zu[] <- 0 ### derivative of a constant
    zl <- x * ftmb[- nrow(ftmb), , drop = FALSE] / prb
    

    ret <- .rowSums(zl - zu, m = nrow(zl), n = ncol(zl)) / 
           .rowSums(x, m = nrow(x), n = ncol(x))
    ret[!is.finite(ret)] <- 0
    return(- ret)
}


.hes <- function(parm, x, mu = 0, rightcensored = FALSE, full = FALSE) 
{

    
    bidx <- seq_len(ncol(x) - 1L)
    delta <- c(0, mu + parm[bidx])
    intercepts <- c(-Inf, parm[- bidx], Inf)
    tmb <- intercepts - matrix(delta, nrow = length(intercepts),  
                                      ncol = ncol(x),
                                      byrow = TRUE)
    Ftmb <- F(tmb)
    if (rightcensored) {
        prb <- 1 - Ftmb[- nrow(Ftmb), , drop = FALSE]
    } else {
        prb <- Ftmb[- 1L, , drop = FALSE] - 
               Ftmb[- nrow(Ftmb), , drop = FALSE]
    } 
    

    
    ftmb <- f(tmb)
    fptmb <- fp(tmb)

    dl <- ftmb[- nrow(ftmb), , drop = FALSE]
    du <- ftmb[- 1, , drop = FALSE]
    if (rightcensored) du[] <- 0
    dpl <- fptmb[- nrow(ftmb), , drop = FALSE]
    dpu <- fptmb[- 1, , drop = FALSE]
    if (rightcensored) dpu[] <- 0
    dlm1 <- dl[,-1L, drop = FALSE]
    dum1 <- du[,-1L, drop = FALSE]
    dplm1 <- dpl[,-1L, drop = FALSE]
    dpum1 <- dpu[,-1L, drop = FALSE]
    prbm1 <- prb[,-1L, drop = FALSE]

    i1 <- length(intercepts) - 1L
    i2 <- 1L
    

    
    Aoffdiag <- - .rowSums(x * du * dl / prb^2, m = nrow(x), n = ncol(x))[-i2]
    Aoffdiag <- Aoffdiag[-length(Aoffdiag)]
    
    
    Adiag <- - .rowSums((x * dpu / prb)[-i1,,drop = FALSE] - 
                        (x * dpl / prb)[-i2,,drop = FALSE] - 
                        ((x * du^2 / prb^2)[-i1,,drop = FALSE] + 
                         (x * dl^2 / prb^2)[-i2,,drop = FALSE] ), 
                        m = nrow(x) - length(i1), n = ncol(x)
                       )
                      
    
    
    xm1 <- x[,-1L,drop = FALSE] 
    X <- ((xm1 * dpum1 / prbm1)[-i1,,drop = FALSE] - 
          (xm1 * dplm1 / prbm1)[-i2,,drop = FALSE] - 
          ((xm1 * dum1^2 / prbm1^2)[-i1,,drop = FALSE] - 
           (xm1 * dum1 * dlm1 / prbm1^2)[-i2,,drop = FALSE] -
           (xm1 * dum1 * dlm1 / prbm1^2)[-i1,,drop = FALSE] +
           (xm1 * dlm1^2 / prbm1^2)[-i2,,drop = FALSE]
          )
         )

    Z <- - .colSums(xm1 * (dpum1 / prbm1 - 
                           dplm1 / prbm1 -
                           (dum1^2 / prbm1^2 - 
                            2 * dum1 * dlm1 / prbm1^2 +
                            dlm1^2 / prbm1^2
                           )
                          ),
                    m = nrow(xm1), n = ncol(xm1)
                    )
    if (length(Z) > 1L) Z <- diag(Z)
    

    if (length(Adiag) > 1L) {
        if (!isFALSE(full)) {
            A <- list(Adiag = Adiag, Aoffdiag = Aoffdiag)
        } else {
            A <- Matrix::bandSparse(length(Adiag), 
                k = 0:1, diagonals = list(Adiag, Aoffdiag), 
                symmetric = TRUE)
        }
    } else {
        if (!isFALSE(full)) {
            A <- list(Adiag = Adiag, Aoffdiag = NULL)
        } else {
            A <- matrix(Adiag)
        }
    }
    return(list(A = A, X = X, Z = Z))
}


.snll <- function(parm, x, mu = 0, rightcensored = FALSE) 
{

    
    C <- vapply(x, NROW, 0L) ### might differ by stratum
    K <- unique(do.call("c", lapply(x, ncol))) ### the same
    B <- length(x)
    sidx <- factor(rep(seq_len(B), times = pmax(0, C - 1L)), 
                   levels = seq_len(B))
    bidx <- seq_len(K - 1L)
    delta <- parm[bidx]
    intercepts <- split(parm[-bidx], sidx)
    

    ret <- 0
    for (b in seq_len(B))
        ret <- ret + .nll(c(delta, intercepts[[b]]), x[[b]], mu = mu,
                          rightcensored = rightcensored)
    return(ret)
}


.snsc <- function(parm, x, mu = 0, rightcensored = FALSE) 
{

    
    C <- vapply(x, NROW, 0L) ### might differ by stratum
    K <- unique(do.call("c", lapply(x, ncol))) ### the same
    B <- length(x)
    sidx <- factor(rep(seq_len(B), times = pmax(0, C - 1L)), 
                   levels = seq_len(B))
    bidx <- seq_len(K - 1L)
    delta <- parm[bidx]
    intercepts <- split(parm[-bidx], sidx)
    

    ret <- numeric(length(bidx))
    for (b in seq_len(B)) {
        nsc <- .nsc(c(delta, intercepts[[b]]), x[[b]], mu = mu,
                    rightcensored = rightcensored)
        ret[bidx] <- ret[bidx] + nsc[bidx]
        ret <- c(ret, nsc[-bidx])
    }
    return(ret)
}


.shes <- function(parm, x, mu = 0, xrc = NULL, full = FALSE, 
                  retMatrix = FALSE) 
{

    
    C <- vapply(x, NROW, 0L) ### might differ by stratum
    K <- unique(do.call("c", lapply(x, ncol))) ### the same
    B <- length(x)
    sidx <- factor(rep(seq_len(B), times = pmax(0, C - 1L)), 
                   levels = seq_len(B))
    bidx <- seq_len(K - 1L)
    delta <- parm[bidx]
    intercepts <- split(parm[-bidx], sidx)
    

    if (!isFALSE(ret <- full)) {
        
        for (b in seq_len(B)) {
            H <- .hes(c(delta, intercepts[[b]]), x[[b]], mu = mu, full = full)
            if (!is.null(xrc)) {
                Hrc <- .hes(c(delta, intercepts[[b]]), xrc[[b]], mu = mu, 
                            rightcensored = TRUE, full = full)
                H$X <- H$X + Hrc$X
                H$A$Adiag <- H$A$Adiag + Hrc$A$Adiag
                H$A$Aoffdiag <- H$A$Aoffdiag + Hrc$A$Aoffdiag
                H$Z <- H$Z + Hrc$Z
            }
            if (b == 1L) {
                Adiag <- H$A$Adiag
                Aoffdiag <- H$A$Aoffdiag
                X <- H$X
                Z <- H$Z
            } else {
                Adiag <- c(Adiag, H$A$Adiag)
                Aoffdiag <- c(Aoffdiag, 0, H$A$Aoffdiag)
                X <- rbind(X, H$X)
                Z <- Z + H$Z
            }
        }

        if (length(Adiag) > 1L) {
            A <- Matrix::bandSparse(length(Adiag),
                                    k = 0:1, diagonals = list(Adiag, Aoffdiag),
                                    symmetric = TRUE)
        } else {
            A <- matrix(Adiag)
        }

        ret <- cbind(Z, t(X))
        ret <- rbind(ret, cbind(X, A))
        if (retMatrix) return(ret)
        return(as.matrix(ret))
        
    }
    ret <- matrix(0, nrow = length(bidx), ncol = length(bidx))
    for (b in seq_len(B)) {
        H <- .hes(c(delta, intercepts[[b]]), x[[b]], mu = mu)
        if (!is.null(xrc)) {
            Hrc <- .hes(c(delta, intercepts[[b]]), xrc[[b]], mu = mu, 
                        rightcensored = TRUE)
            H$X <- H$X + Hrc$X
            H$A <- H$A + Hrc$A
            H$Z <- H$Z + Hrc$Z
        }
        sAH <- tryCatch(Matrix::solve(H$A, H$X), error = function(e) NULL)
        if (is.null(sAH))
            stop(gettextf("error computing the Hessian in %s",
                          "free1way"),
                 domain = NA)
        ret <- ret + (H$Z - crossprod(H$X, sAH))
    }
    as.matrix(ret)
}


.snsr <- function(parm, x, mu = 0, rightcensored = FALSE) 
{

    
    C <- vapply(x, NROW, 0L) ### might differ by stratum
    K <- unique(do.call("c", lapply(x, ncol))) ### the same
    B <- length(x)
    sidx <- factor(rep(seq_len(B), times = pmax(0, C - 1L)), 
                   levels = seq_len(B))
    bidx <- seq_len(K - 1L)
    delta <- parm[bidx]
    intercepts <- split(parm[-bidx], sidx)
    

    ret <- c()
    for (b in seq_len(B)) {
        idx <- attr(x[[b]], "idx")
        ### idx == 1L means zero residual, see definition of idx
        sr <- c(0, .nsr(c(delta, intercepts[[b]]), x[[b]], mu = mu,
                        rightcensored = rightcensored))
        ret <- c(ret, sr[idx])
    }
    return(ret)
}


.free1wayML <- function(x, link, mu = 0, start = NULL, fix = NULL, 
                        residuals = TRUE, score = TRUE, hessian = TRUE, 
                        MPL_Jeffreys = FALSE,
                        ### use nlminb for small sample sizes
                        dooptim = c(".NewtonRaphson", "nlminb")[1 + (sum(x) < 20)],                         
                        control = list(
                            "nlminb" = list(trace = trace, iter.max = 200,
                                            eval.max = 200, rel.tol = 1e-10,
                                            abs.tol = 1e-20, xf.tol = 1e-16),
                            ".NewtonRaphson" = list(iter.max = 200, trace = trace, 
                                             objtol = 5e-4, 
                                             gradtol = 1e-5 * sum(x) / 1000, 
                                             paramtol = 1e-5, minstepsize = 1e-2, 
                                             tolsolve = .Machine$double.eps)
                        )[dooptim],
                        trace = FALSE, 
                        tol = sqrt(.Machine$double.eps), ...) 
{

    ### convert to three-way table
    xt <- x
    if (!is.table(x))
        stop(gettextf("invalid argument '%s'", "x"), domain = NA) # 'y' in free1way ...
    dx <- dim(x)
    dn <- dimnames(x)
    if (length(dx) == 2L) {
        x <- as.table(array(c(x), dim = dx <- c(dx, 1L)))
        dimnames(x) <- dn <- c(dn, list(A = "A"))
    }

    ### short-cuts for link functions
    F <- function(q) .p(link, q = q)
    Q <- function(p) .q(link, p = p)
    f <- function(q) .d(link, x = q)
    fp <- function(q) .dd(link, x = q)

    if(!suppressPackageStartupMessages(requireNamespace("Matrix")))
        stop(gettextf("%s needs package 'Matrix' correctly installed",
                      ".free1wayML"),
                 domain = NA)

    

    dx <- dim(x)
    if (length(dx) == 1L)
        stop("incorrect dimensions")
    if (length(dx) == 2L)
        x <- as.table(array(x, dim = c(dx, 1)))
    dx <- dim(x)
    if (length(dx) < 3L)
        stop("incorrect dimensions")
    C <- dim(x)[1L]
    K <- dim(x)[2L]
    B <- dim(x)[3L]
    if (C < 2L)
        stop("at least two response categories required")
    if (K < 2L)
        stop("at least two groups required")
    xrc <- NULL
    if (length(dx) == 4L) {
        if (dx[4] == 2L) {
            xrc <- array(x[,,,"FALSE", drop = TRUE], dim = dx[1:3])
            x <- array(x[,,,"TRUE", drop = TRUE], dim = dx[1:3])
        } else {
            stop(gettextf("%s currently only allows independent right-censoring",
                          "free1way"),
                 domain = NA)
        }
    }


    xlist <- xrclist <- vector(mode = "list", length = B)

    for (b in seq_len(B)) {
        xb <- matrix(x[,,b, drop = TRUE], ncol = K)
        xw <- rowSums(abs(xb)) > 0
        if (sum(xw) > 1L) {
            ### do not remove last parameter if there are corresponding
            ### right-censored observations
            wm <- which(xw)[sum(xw)]
            if (!is.null(xrc) && any(xrc[wm:dx[1],,b,drop = TRUE] > 0))
                xw[length(xw)] <- TRUE
            xlist[[b]] <- xb[xw,,drop = FALSE]
            Cidx <- rep.int(1L, times = C)
            Cidx[xw] <- Cidx[xw] + seq_len(sum(xw))
            attr(xlist[[b]], "idx") <- Cidx
            if (!is.null(xrc)) {
                ### count right-censored observations between distinct event
                ### times
                cs <- apply(xrc[,,b,drop = TRUE] * (!xw), 2, function(x) 
                    diff(c(0, cumsum(x)[xw])))
                xrclist[[b]] <- matrix(xrc[xw,,b,drop = TRUE], ncol = K) + cs
                idx <- seq_len(C)[xw]
                idx <- rep(seq_len(sum(xw)), times = c(idx[1], diff(idx)))
                Cidx <- rep.int(1L, times = C)
                Cidx[seq_along(idx)] <- Cidx[seq_along(idx)] + idx
                attr(xrclist[[b]], "idx") <- Cidx
            }
        }
    }
    ### remove empty blocks
    strata <- !vapply(xlist, is.null, NA)
    xlist <- xlist[strata]
    xrclist <- xrclist[strata]
    
    
    ## allow specification of start = delta and fix = 1:K
    ## for evaluating the likelihood at given delta parameters
    ## without having to specify all intercept parameters
    if (is.null(start))
        start <- rep.int(0, K - 1L)
    NS <- length(start) == (K - 1L)
    lwr <- rep(-Inf, times = K - 1L)
    for (b in seq_len(length(xlist))) {
        bC <- nrow(xlist[[b]]) - 1L
        lwr <- c(lwr, -Inf, rep.int(0, times = bC - 1L))
        if (NS) {
            ecdf0 <- cumsum(rowSums(xlist[[b]]))
            ### ensure that 0 < ecdf0 < 1 such that quantiles exist
            ecdf0 <- pmax(1, ecdf0[-length(ecdf0)]) / (max(ecdf0) + 1)
            Qecdf <- Q(ecdf0)
            start <- c(start, Qecdf)
        }
    }
    
    
    .nll <- function(parm, x, mu = 0, rightcensored = FALSE) 
    {

        
        bidx <- seq_len(ncol(x) - 1L)
        delta <- c(0, mu + parm[bidx])
        intercepts <- c(-Inf, parm[- bidx], Inf)
        tmb <- intercepts - matrix(delta, nrow = length(intercepts),  
                                          ncol = ncol(x),
                                          byrow = TRUE)
        Ftmb <- F(tmb)
        if (rightcensored) {
            prb <- 1 - Ftmb[- nrow(Ftmb), , drop = FALSE]
        } else {
            prb <- Ftmb[- 1L, , drop = FALSE] - 
                   Ftmb[- nrow(Ftmb), , drop = FALSE]
        } 
        
        if (any(prb < .Machine$double.eps^10)) 
            return(Inf)
        return(- sum(x * log(prb)))
    }
    
    
    .nsc <- function(parm, x, mu = 0, rightcensored = FALSE) 
    {

        
        bidx <- seq_len(ncol(x) - 1L)
        delta <- c(0, mu + parm[bidx])
        intercepts <- c(-Inf, parm[- bidx], Inf)
        tmb <- intercepts - matrix(delta, nrow = length(intercepts),  
                                          ncol = ncol(x),
                                          byrow = TRUE)
        Ftmb <- F(tmb)
        if (rightcensored) {
            prb <- 1 - Ftmb[- nrow(Ftmb), , drop = FALSE]
        } else {
            prb <- Ftmb[- 1L, , drop = FALSE] - 
                   Ftmb[- nrow(Ftmb), , drop = FALSE]
        } 
        

        
        ftmb <- f(tmb)
        zu <- x * ftmb[- 1, , drop = FALSE] / prb
        if (rightcensored) zu[] <- 0 ### derivative of a constant
        zl <- x * ftmb[- nrow(ftmb), , drop = FALSE] / prb
        

        ret <- numeric(length(parm))
        ret[bidx] <- .colSums(zl, m = nrow(zl), n = ncol(zl))[-1L] -
                     .colSums(zu[-nrow(zu),,drop = FALSE], 
                              m = nrow(zu) - 1L, n = ncol(zu))[-1L]
        ret[- bidx] <- .rowSums(zu[-nrow(zu),,drop = FALSE] - 
                                zl[-1,,drop = FALSE], 
                                m = nrow(zu) - 1L, n = ncol(zu))
        return(- ret)
    }
    
    
    .nsr <- function(parm, x, mu = 0, rightcensored = FALSE) 
    {

        
        bidx <- seq_len(ncol(x) - 1L)
        delta <- c(0, mu + parm[bidx])
        intercepts <- c(-Inf, parm[- bidx], Inf)
        tmb <- intercepts - matrix(delta, nrow = length(intercepts),  
                                          ncol = ncol(x),
                                          byrow = TRUE)
        Ftmb <- F(tmb)
        if (rightcensored) {
            prb <- 1 - Ftmb[- nrow(Ftmb), , drop = FALSE]
        } else {
            prb <- Ftmb[- 1L, , drop = FALSE] - 
                   Ftmb[- nrow(Ftmb), , drop = FALSE]
        } 
        

        
        ftmb <- f(tmb)
        zu <- x * ftmb[- 1, , drop = FALSE] / prb
        if (rightcensored) zu[] <- 0 ### derivative of a constant
        zl <- x * ftmb[- nrow(ftmb), , drop = FALSE] / prb
        

        ret <- .rowSums(zl - zu, m = nrow(zl), n = ncol(zl)) / 
               .rowSums(x, m = nrow(x), n = ncol(x))
        ret[!is.finite(ret)] <- 0
        return(- ret)
    }
    
    
    .hes <- function(parm, x, mu = 0, rightcensored = FALSE, full = FALSE) 
    {

        
        bidx <- seq_len(ncol(x) - 1L)
        delta <- c(0, mu + parm[bidx])
        intercepts <- c(-Inf, parm[- bidx], Inf)
        tmb <- intercepts - matrix(delta, nrow = length(intercepts),  
                                          ncol = ncol(x),
                                          byrow = TRUE)
        Ftmb <- F(tmb)
        if (rightcensored) {
            prb <- 1 - Ftmb[- nrow(Ftmb), , drop = FALSE]
        } else {
            prb <- Ftmb[- 1L, , drop = FALSE] - 
                   Ftmb[- nrow(Ftmb), , drop = FALSE]
        } 
        

        
        ftmb <- f(tmb)
        fptmb <- fp(tmb)

        dl <- ftmb[- nrow(ftmb), , drop = FALSE]
        du <- ftmb[- 1, , drop = FALSE]
        if (rightcensored) du[] <- 0
        dpl <- fptmb[- nrow(ftmb), , drop = FALSE]
        dpu <- fptmb[- 1, , drop = FALSE]
        if (rightcensored) dpu[] <- 0
        dlm1 <- dl[,-1L, drop = FALSE]
        dum1 <- du[,-1L, drop = FALSE]
        dplm1 <- dpl[,-1L, drop = FALSE]
        dpum1 <- dpu[,-1L, drop = FALSE]
        prbm1 <- prb[,-1L, drop = FALSE]

        i1 <- length(intercepts) - 1L
        i2 <- 1L
        

        
        Aoffdiag <- - .rowSums(x * du * dl / prb^2, m = nrow(x), n = ncol(x))[-i2]
        Aoffdiag <- Aoffdiag[-length(Aoffdiag)]
        
        
        Adiag <- - .rowSums((x * dpu / prb)[-i1,,drop = FALSE] - 
                            (x * dpl / prb)[-i2,,drop = FALSE] - 
                            ((x * du^2 / prb^2)[-i1,,drop = FALSE] + 
                             (x * dl^2 / prb^2)[-i2,,drop = FALSE] ), 
                            m = nrow(x) - length(i1), n = ncol(x)
                           )
                          
        
        
        xm1 <- x[,-1L,drop = FALSE] 
        X <- ((xm1 * dpum1 / prbm1)[-i1,,drop = FALSE] - 
              (xm1 * dplm1 / prbm1)[-i2,,drop = FALSE] - 
              ((xm1 * dum1^2 / prbm1^2)[-i1,,drop = FALSE] - 
               (xm1 * dum1 * dlm1 / prbm1^2)[-i2,,drop = FALSE] -
               (xm1 * dum1 * dlm1 / prbm1^2)[-i1,,drop = FALSE] +
               (xm1 * dlm1^2 / prbm1^2)[-i2,,drop = FALSE]
              )
             )

        Z <- - .colSums(xm1 * (dpum1 / prbm1 - 
                               dplm1 / prbm1 -
                               (dum1^2 / prbm1^2 - 
                                2 * dum1 * dlm1 / prbm1^2 +
                                dlm1^2 / prbm1^2
                               )
                              ),
                        m = nrow(xm1), n = ncol(xm1)
                        )
        if (length(Z) > 1L) Z <- diag(Z)
        

        if (length(Adiag) > 1L) {
            if (!isFALSE(full)) {
                A <- list(Adiag = Adiag, Aoffdiag = Aoffdiag)
            } else {
                A <- Matrix::bandSparse(length(Adiag), 
                    k = 0:1, diagonals = list(Adiag, Aoffdiag), 
                    symmetric = TRUE)
            }
        } else {
            if (!isFALSE(full)) {
                A <- list(Adiag = Adiag, Aoffdiag = NULL)
            } else {
                A <- matrix(Adiag)
            }
        }
        return(list(A = A, X = X, Z = Z))
    }
    
    
    .snll <- function(parm, x, mu = 0, rightcensored = FALSE) 
    {

        
        C <- vapply(x, NROW, 0L) ### might differ by stratum
        K <- unique(do.call("c", lapply(x, ncol))) ### the same
        B <- length(x)
        sidx <- factor(rep(seq_len(B), times = pmax(0, C - 1L)), 
                       levels = seq_len(B))
        bidx <- seq_len(K - 1L)
        delta <- parm[bidx]
        intercepts <- split(parm[-bidx], sidx)
        

        ret <- 0
        for (b in seq_len(B))
            ret <- ret + .nll(c(delta, intercepts[[b]]), x[[b]], mu = mu,
                              rightcensored = rightcensored)
        return(ret)
    }
    
    
    .snsc <- function(parm, x, mu = 0, rightcensored = FALSE) 
    {

        
        C <- vapply(x, NROW, 0L) ### might differ by stratum
        K <- unique(do.call("c", lapply(x, ncol))) ### the same
        B <- length(x)
        sidx <- factor(rep(seq_len(B), times = pmax(0, C - 1L)), 
                       levels = seq_len(B))
        bidx <- seq_len(K - 1L)
        delta <- parm[bidx]
        intercepts <- split(parm[-bidx], sidx)
        

        ret <- numeric(length(bidx))
        for (b in seq_len(B)) {
            nsc <- .nsc(c(delta, intercepts[[b]]), x[[b]], mu = mu,
                        rightcensored = rightcensored)
            ret[bidx] <- ret[bidx] + nsc[bidx]
            ret <- c(ret, nsc[-bidx])
        }
        return(ret)
    }
    
    
    .shes <- function(parm, x, mu = 0, xrc = NULL, full = FALSE, 
                      retMatrix = FALSE) 
    {

        
        C <- vapply(x, NROW, 0L) ### might differ by stratum
        K <- unique(do.call("c", lapply(x, ncol))) ### the same
        B <- length(x)
        sidx <- factor(rep(seq_len(B), times = pmax(0, C - 1L)), 
                       levels = seq_len(B))
        bidx <- seq_len(K - 1L)
        delta <- parm[bidx]
        intercepts <- split(parm[-bidx], sidx)
        

        if (!isFALSE(ret <- full)) {
            
            for (b in seq_len(B)) {
                H <- .hes(c(delta, intercepts[[b]]), x[[b]], mu = mu, full = full)
                if (!is.null(xrc)) {
                    Hrc <- .hes(c(delta, intercepts[[b]]), xrc[[b]], mu = mu, 
                                rightcensored = TRUE, full = full)
                    H$X <- H$X + Hrc$X
                    H$A$Adiag <- H$A$Adiag + Hrc$A$Adiag
                    H$A$Aoffdiag <- H$A$Aoffdiag + Hrc$A$Aoffdiag
                    H$Z <- H$Z + Hrc$Z
                }
                if (b == 1L) {
                    Adiag <- H$A$Adiag
                    Aoffdiag <- H$A$Aoffdiag
                    X <- H$X
                    Z <- H$Z
                } else {
                    Adiag <- c(Adiag, H$A$Adiag)
                    Aoffdiag <- c(Aoffdiag, 0, H$A$Aoffdiag)
                    X <- rbind(X, H$X)
                    Z <- Z + H$Z
                }
            }

            if (length(Adiag) > 1L) {
                A <- Matrix::bandSparse(length(Adiag),
                                        k = 0:1, diagonals = list(Adiag, Aoffdiag),
                                        symmetric = TRUE)
            } else {
                A <- matrix(Adiag)
            }

            ret <- cbind(Z, t(X))
            ret <- rbind(ret, cbind(X, A))
            if (retMatrix) return(ret)
            return(as.matrix(ret))
            
        }
        ret <- matrix(0, nrow = length(bidx), ncol = length(bidx))
        for (b in seq_len(B)) {
            H <- .hes(c(delta, intercepts[[b]]), x[[b]], mu = mu)
            if (!is.null(xrc)) {
                Hrc <- .hes(c(delta, intercepts[[b]]), xrc[[b]], mu = mu, 
                            rightcensored = TRUE)
                H$X <- H$X + Hrc$X
                H$A <- H$A + Hrc$A
                H$Z <- H$Z + Hrc$Z
            }
            sAH <- tryCatch(Matrix::solve(H$A, H$X), error = function(e) NULL)
            if (is.null(sAH))
                stop(gettextf("error computing the Hessian in %s",
                              "free1way"),
                     domain = NA)
            ret <- ret + (H$Z - crossprod(H$X, sAH))
        }
        as.matrix(ret)
    }
    
    
    .snsr <- function(parm, x, mu = 0, rightcensored = FALSE) 
    {

        
        C <- vapply(x, NROW, 0L) ### might differ by stratum
        K <- unique(do.call("c", lapply(x, ncol))) ### the same
        B <- length(x)
        sidx <- factor(rep(seq_len(B), times = pmax(0, C - 1L)), 
                       levels = seq_len(B))
        bidx <- seq_len(K - 1L)
        delta <- parm[bidx]
        intercepts <- split(parm[-bidx], sidx)
        

        ret <- c()
        for (b in seq_len(B)) {
            idx <- attr(x[[b]], "idx")
            ### idx == 1L means zero residual, see definition of idx
            sr <- c(0, .nsr(c(delta, intercepts[[b]]), x[[b]], mu = mu,
                            rightcensored = rightcensored))
            ret <- c(ret, sr[idx])
        }
        return(ret)
    }
    
    

    fn <- function(par) 
    {
        ret <- .snll(par, x = xlist, mu = mu)
        if (!is.null(xrc))
            ret <- ret + .snll(par, x = xrclist, mu = mu, 
                               rightcensored = TRUE)
        return(ret)
    }
    gr <- function(par) 
    {
        ret <- .snsc(par, x = xlist, mu = mu)
        if (!is.null(xrc))
            ret <- ret + .snsc(par, x = xrclist, mu = mu, 
                               rightcensored = TRUE)
        return(ret)
    }

    ### allocate memory for hessian
    Hess <- Matrix::Matrix(0, nrow = length(start), ncol = length(start))

    he <- function(par) 
    {
        if (!is.null(xrc)) {
            ret <- .shes(par, x = xlist, mu = mu, xrc = xrclist, full = Hess, 
                         retMatrix = names(control)[1L] == ".NewtonRaphson")
        } else {
            ret <- .shes(par, x = xlist, mu = mu, full = Hess, 
                         retMatrix = names(control)[1L] == ".NewtonRaphson")
        }
        return(ret)
    }
    

    .profile <- function(start, fix = seq_len(K - 1)) 
    {
        if (!all(fix %in% seq_len(K - 1)))
            stop(gettextf("invalid argument '%s'", "fix"), domain = NA)
        delta <- start[fix]
        opargs <- list(start = start[-fix], 
                         objective = function(par) {
                             p <- numeric(length(par) + length(fix))
                             p[fix] <- delta
                             p[-fix] <- par
                             fn(p)
                         },
                         gradient = function(par) {
                             p <- numeric(length(par) + length(fix))
                             p[fix] <- delta
                             p[-fix] <- par
                             gr(p)[-fix]
                         },
                         hessian = function(par) {
                             p <- numeric(length(par) + length(fix))
                             p[fix] <- delta
                             p[-fix] <- par
                             he(p)[-fix, -fix, drop = FALSE]
                         })
        opargs$control <- control[[1L]]
        MPL_Jeffreys <- FALSE ### turn off Jeffreys penalisation in .profile

        
        maxit <- control[[1L]]$iter.max
        while(maxit < 10001) {
           ret <- do.call(names(control)[[1L]], opargs)
           maxit <- 5 * maxit
           if (ret$convergence > 0) {
               opargs$control$eval.max <- maxit
               opargs$control$iter.max <- maxit
               opargs$start <- ret$par
           } else {
               break()
           }
        }

        if (isTRUE(MPL_Jeffreys)) {
            
            .pll_Jeffreys <- function(cf, start) 
            {
                fix <- seq_along(cf)
                start[fix] <- cf
                ### compute profile likelihood w/o warnings
                ret <- suppressWarnings(.profile(start, fix = fix))
                Hfull <- he(ret$par)
                Hfix <- as.matrix(solve(solve(Hfull)[fix, fix]))
                return(ret$value - 
                       .5 * determinant(Hfix, logarithm = TRUE)$modulus)
            }
            if (K == 2) {
                MLcf <- ret$par[seq_len(K - 1)]
                Fret <- optim(MLcf, fn = .pll_Jeffreys, start = ret$par,
                              method = "Brent", lower = MLcf - 5, 
                              upper = MLcf + 5)
            } else {
                ### Nelder-Mead
                Fret <- optim(ret$par[seq_len(K - 1)], fn = .pll_Jeffreys, 
                              start = ret$par)
            }
            if (Fret$convergence == 0) {
                start <- ret$par
                start[seq_len(K - 1)] <- Fret$par
                ret <- .profile(start, fix = seq_len(K - 1))
                ret$objective <- ret$value
            }
            
        } else {
            if (ret$convergence > 0) {
                if (is.na(MPL_Jeffreys)) { ### only after failure
                    warning(gettextf("Jeffreys penalisation was applied in %s because initial optimisation failed with:",
                                     "free1way"),
                            "\n  ", ret$message, domain = NA)
                    MPL_Jeffreys <- TRUE
                    
                    .pll_Jeffreys <- function(cf, start) 
                    {
                        fix <- seq_along(cf)
                        start[fix] <- cf
                        ### compute profile likelihood w/o warnings
                        ret <- suppressWarnings(.profile(start, fix = fix))
                        Hfull <- he(ret$par)
                        Hfix <- as.matrix(solve(solve(Hfull)[fix, fix]))
                        return(ret$value - 
                               .5 * determinant(Hfix, logarithm = TRUE)$modulus)
                    }
                    if (K == 2) {
                        MLcf <- ret$par[seq_len(K - 1)]
                        Fret <- optim(MLcf, fn = .pll_Jeffreys, start = ret$par,
                                      method = "Brent", lower = MLcf - 5, 
                                      upper = MLcf + 5)
                    } else {
                        ### Nelder-Mead
                        Fret <- optim(ret$par[seq_len(K - 1)], fn = .pll_Jeffreys, 
                                      start = ret$par)
                    }
                    if (Fret$convergence == 0) {
                        start <- ret$par
                        start[seq_len(K - 1)] <- Fret$par
                        ret <- .profile(start, fix = seq_len(K - 1))
                        ret$objective <- ret$value
                    }
                    
                }
           }
        }
        if (ret$convergence > 0)
            warning(gettextf("unsuccessful optimisation in %s", "free1way"),
                    ": ", ret$message, domain = NA)

        ret$MPL_Jeffreys <- MPL_Jeffreys
        ret$value <- ret$objective
        ret$objective <- NULL
        

        p <- numeric(length(start))
        p[fix] <- delta
        p[-fix] <- ret$par
        ret$par <- p
        ret
    }
    
    
    if (!length(fix)) {
        opargs <- list(start = start, 
                       objective = fn, 
                       gradient = gr,
                       hessian = he)
        opargs$control <- control[[1L]]
        
        maxit <- control[[1L]]$iter.max
        while(maxit < 10001) {
           ret <- do.call(names(control)[[1L]], opargs)
           maxit <- 5 * maxit
           if (ret$convergence > 0) {
               opargs$control$eval.max <- maxit
               opargs$control$iter.max <- maxit
               opargs$start <- ret$par
           } else {
               break()
           }
        }

        if (isTRUE(MPL_Jeffreys)) {
            
            .pll_Jeffreys <- function(cf, start) 
            {
                fix <- seq_along(cf)
                start[fix] <- cf
                ### compute profile likelihood w/o warnings
                ret <- suppressWarnings(.profile(start, fix = fix))
                Hfull <- he(ret$par)
                Hfix <- as.matrix(solve(solve(Hfull)[fix, fix]))
                return(ret$value - 
                       .5 * determinant(Hfix, logarithm = TRUE)$modulus)
            }
            if (K == 2) {
                MLcf <- ret$par[seq_len(K - 1)]
                Fret <- optim(MLcf, fn = .pll_Jeffreys, start = ret$par,
                              method = "Brent", lower = MLcf - 5, 
                              upper = MLcf + 5)
            } else {
                ### Nelder-Mead
                Fret <- optim(ret$par[seq_len(K - 1)], fn = .pll_Jeffreys, 
                              start = ret$par)
            }
            if (Fret$convergence == 0) {
                start <- ret$par
                start[seq_len(K - 1)] <- Fret$par
                ret <- .profile(start, fix = seq_len(K - 1))
                ret$objective <- ret$value
            }
            
        } else {
            if (ret$convergence > 0) {
                if (is.na(MPL_Jeffreys)) { ### only after failure
                    warning(gettextf("Jeffreys penalisation was applied in %s because initial optimisation failed with:",
                                     "free1way"),
                            "\n  ", ret$message, domain = NA)
                    MPL_Jeffreys <- TRUE
                    
                    .pll_Jeffreys <- function(cf, start) 
                    {
                        fix <- seq_along(cf)
                        start[fix] <- cf
                        ### compute profile likelihood w/o warnings
                        ret <- suppressWarnings(.profile(start, fix = fix))
                        Hfull <- he(ret$par)
                        Hfix <- as.matrix(solve(solve(Hfull)[fix, fix]))
                        return(ret$value - 
                               .5 * determinant(Hfix, logarithm = TRUE)$modulus)
                    }
                    if (K == 2) {
                        MLcf <- ret$par[seq_len(K - 1)]
                        Fret <- optim(MLcf, fn = .pll_Jeffreys, start = ret$par,
                                      method = "Brent", lower = MLcf - 5, 
                                      upper = MLcf + 5)
                    } else {
                        ### Nelder-Mead
                        Fret <- optim(ret$par[seq_len(K - 1)], fn = .pll_Jeffreys, 
                                      start = ret$par)
                    }
                    if (Fret$convergence == 0) {
                        start <- ret$par
                        start[seq_len(K - 1)] <- Fret$par
                        ret <- .profile(start, fix = seq_len(K - 1))
                        ret$objective <- ret$value
                    }
                    
                }
           }
        }
        if (ret$convergence > 0)
            warning(gettextf("unsuccessful optimisation in %s", "free1way"),
                    ": ", ret$message, domain = NA)

        ret$MPL_Jeffreys <- MPL_Jeffreys
        ret$value <- ret$objective
        ret$objective <- NULL
        
    } else if (length(fix) == length(start)) {
        ret <- list(par = start, 
                    value = fn(start))
    } else {
        ret <- .profile(start, fix = fix)
    }
     
    
    if (is.null(fix) || (length(fix) == length(start)))
        parm <- seq_len(K - 1)
    else 
        parm <- fix
    if (any(parm >= K)) return(ret)

    ret$coefficients <- ret$par[parm]
    dn2 <- dimnames(xt)[2L]
    names(ret$coefficients) <- cnames <- paste0(names(dn2), dn2[[1L]][1L + parm])

    par <- ret$par
    intercepts <- function(parm, x) 
    {

        
        C <- vapply(x, NROW, 0L) ### might differ by stratum
        K <- unique(do.call("c", lapply(x, ncol))) ### the same
        B <- length(x)
        sidx <- factor(rep(seq_len(B), times = pmax(0, C - 1L)), 
                       levels = seq_len(B))
        bidx <- seq_len(K - 1L)
        delta <- parm[bidx]
        intercepts <- split(parm[-bidx], sidx)
        

        return(intercepts)
    }
    ret$intercepts <- intercepts(par, x = xlist)

    if (score) {
        ret$negscore <- .snsc(par, x = xlist, mu = mu)[parm]
        if (!is.null(xrc))
            ret$negscore <- ret$negscore + .snsc(par, x = xrclist, mu = mu, 
                                                 rightcensored = TRUE)[parm]
    }
    if (hessian) {
        if (!is.null(xrc)) {
            ret$hessian <- .shes(par, x = xlist, mu = mu, xrc = xrclist)
        } else {
            ret$hessian <- .shes(par, x = xlist, mu = mu)
        }
        ret$vcov <- solve(ret$hessian)
        if (length(parm) != nrow(ret$hessian))
           ret$hessian <- solve(ret$vcov <- ret$vcov[parm, parm, drop = FALSE])
        rownames(ret$vcov) <- colnames(ret$vcov) <- rownames(ret$hessian) <-
            colnames(ret$hessian) <-  cnames
    }
    if (residuals) {
        ret$negresiduals <- .snsr(par, x = xlist, mu = mu)
        if (!is.null(xrc)) {
            rcr <- .snsr(par, x = xrclist, mu = mu, rightcensored = TRUE)
            ret$negresiduals <- c(rbind(matrix(ret$negresiduals, nrow = C),
                                        matrix(rcr, nrow = C)))
         }
    }
    ret$profile <- function(start, fix)
        .free1wayML(xt, link = link, mu = mu, start = start, fix = fix, tol = tol, 
                    ...) 
    ret$table <- xt

    ret$strata <- strata
    ret$mu <- mu
    if (length(ret$mu) == 1) {
        names(ret$mu) <- link$parm
    } else {
        names(ret$mu) <- c(paste(link$parm, cnames[1L], sep = ":"), cnames[-1L])
    }
    

    class(ret) <- "free1wayML"
    ret
}


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


### distribution function
.p <- function(link, q, ...)
    link$linkinv(q = q, ...)

### quantile function
.q <- function(link, p, ...)
    link$link(p = p, ...)

### density function
.d <- function(link, x, ...)
    link$dlinkinv(x = x, ...)

### derivative of density function
.dd <- function(link, x, ...)
    link$ddlinkinv(x = x, ...)

### 2nd derivative of density function
.ddd <- function(link, x, ...)
    link$dddlinkinv(x = x, ...)

### ratio of derivative of density to
### density function
.dd2d <- function(link, x, ...)
    link$dd2dlinkinv(x = x, ...)

### constructor
linkfun <- function(name,	### nickname
                    alias,	### char 
                    model, 	### char, semiparametric model name
                    parm, 	### char, parameter name
                    link,      	### quantile function
                    linkinv,   	### distribution function
                    dlinkinv,  	### density function
                    ddlinkinv, 	### derivative of density function
                    ...) 
{

    ret <- list(name = name, 
                alias = alias,
                model = model,
                parm = parm,
                link = link,
                linkinv = linkinv,
                dlinkinv = dlinkinv,
                ddlinkinv = ddlinkinv)
    if (is.null(ret$dd2d)) 
        ret$dd2d <- function(x) 
            ret$ddlinkinv(x) / ret$dlinkinv(x)
    ret <- c(ret, list(...))
    class(ret) <- "linkfun"
    ret
}


logit <- function()
    linkfun(name = "Logit", 
            alias = c("Wilcoxon", "Kruskal-Wallis"),
            model = "proportional odds", 
            parm = "log-odds ratio",
            link = qlogis,
            linkinv = plogis,
            dlinkinv = dlogis,
            ddlinkinv = function(x) {
                p <- plogis(x)
                p * (1 - p)^2 - p^2 * (1 - p)
            },
            dddlinkinv = function(x) {
                ex <- exp(x)
                ifelse(is.finite(x), (ex - 4 * ex^2 + ex^3) / (1 + ex)^4, 0.0)
            },
            dd2d = function(x) {
                ex <- exp(x)
                (1 - ex) / (1 + ex)
            },
            parm2PI = function(x) {
               OR <- exp(x)
               ret <- OR * (OR - 1 - x)/(OR - 1)^2
               ret[abs(x) < .Machine$double.eps] <- 0.5
               return(ret)
            },
            PI2parm = function(p) {
               f <- function(x, PI)
                   x + (exp(-x) * (PI + 
                                   exp(2 * x) * (PI - 1) + 
                                   exp(x) * (1 - 2 * PI)))
               ret <- vapply(p, function(p) 
                   uniroot(f, PI = p, interval = 50 * c(-1, 1))$root, 0)
               return(ret)
            },
            parm2OVL = function(x) 2 * plogis(-abs(x / 2))
    )


probit <- function()
    linkfun(name = "Probit",
            alias = "van der Waerden normal scores",
            model = "latent normal shift", 
            parm = "generalised Cohen's d",
            link = qnorm,
            linkinv = pnorm,
            dlinkinv = dnorm,
            ddlinkinv = function(x) 
                ifelse(is.finite(x), -dnorm(x = x) * x, 0.0), 
            dddlinkinv = function(x) 
                ifelse(is.finite(x), dnorm(x = x) * (x^2 - 1), 0.0),
            dd2d = function(x) -x,
            parm2PI = function(x) pnorm(x, sd = sqrt(2)),
            PI2parm = function(p) qnorm(p, sd = sqrt(2)),
            parm2OVL = function(x) 2 * pnorm(-abs(x / 2))
    )


cloglog <- function()
    linkfun(name = "Complementary Log-log",
            alias = "Savage",
            model = "proportional hazards", 
            parm = "log-hazard ratio",
            link = function(p, log.p = FALSE) {
                if (log.p) p <- exp(p)
                log(-log1p(- p))
            },
            linkinv = function(q, lower.tail = TRUE, log.p = FALSE) {
                ### p = 1 - exp(-exp(q))
                ret <- exp(-exp(q))
                if (log.p) {
                    if (lower.tail)
                        return(log1p(-ret))
                    return(-exp(q))
                }
                if (lower.tail)
                    return(-expm1(-exp(q)))
                return(ret)
            },
            dlinkinv = function(x) 
                ifelse(is.finite(x), exp(x - exp(x)), 0.0),
            ddlinkinv = function(x) {
                ex <- exp(x)
                ifelse(is.finite(x), (ex - ex^2) / exp(ex), 0.0)
            },
            dddlinkinv = function(x) {
                ex <- exp(x)
                ifelse(is.finite(x), (ex - 3*ex^2 + ex^3) / exp(ex), 0.0)
            },
            dd2d = function(x)
               -expm1(x),
            parm2PI = plogis,
            PI2parm = qlogis,
            parm2OVL = function(x) {
                x <- abs(x)
                ret <- exp(x / (exp(-x) - 1)) - exp(-x / (exp(x) - 1)) + 1 
                ret[abs(x) < .Machine$double.eps] <- 1
                x[] <- ret
                return(x)
            }
    )


loglog <- function()
    linkfun(name = "Log-log",
            alias = "Lehmann", 
            model = "Lehmann", 
            parm = "log-reverse time hazard ratio",
            link = function(p, log.p = FALSE) {
                if (!log.p) p <- log(p)
                -log(-p)
            },
            linkinv = function(q, lower.tail = TRUE, log.p = FALSE) {
                ### p = exp(-exp(-q))
                if (log.p) {
                    if (lower.tail)
                        return(-exp(-q))
                    return(log1p(-exp(-exp(-q))))
                }
                if (lower.tail)
                    return(exp(-exp(-q)))
                -expm1(-exp(-q))
            },
            dlinkinv = function(x) 
                ifelse(is.finite(x), exp(- x - exp(-x)), 0.0),
            ddlinkinv = function(x) {
               ex <- exp(-x)
               ifelse(is.finite(x), exp(-ex - x) * (ex - 1.0), 0.0)
            },
            dddlinkinv = function(x) {
               ex <- exp(-x)
               ifelse(is.finite(x), exp(-x - ex) * (ex - 1)^2 - 
                                    exp(-ex - 2 * x), 
                                    0.0)
            },
            dd2d = function(x) 
                expm1(-x),
            parm2PI = plogis,
            PI2parm = qlogis,
            parm2OVL = function(x) {
                x <- abs(x)
                rt <- exp(-x / (exp(x) - 1))
                ret <- rt^exp(x) + 1 - rt
                ret[abs(x) < .Machine$double.eps] <- 1
                x[] <- ret
                return(x)
            }
    )


### adopted from rms:::lrm.fit
.NewtonRaphson <- function(start, objective, gradient, hessian, 
                           control = list(iter.max = 150, trace = trace, 
                                          objtol = 5e-4, gradtol = 1e-5, 
                                          paramtol = 1e-5, minstepsize = 1e-2, 
                                          tolsolve = .Machine$double.eps),
                           trace = FALSE)
{

    theta  <- start # Initialize the parameter vector
    oldobj <- Inf
    objthe <- objective(theta)
    if (!is.finite(objthe)) {
        msg <- "Infeasible starting values"
        return(list(par = theta, objective = objthe, convergence = 1, 
                    message = msg)) 
    }

    if(!suppressPackageStartupMessages(requireNamespace("Matrix")))
        stop(gettextf("%s needs package 'Matrix' correctly installed",
                      ".NewtonRaphson"),
                 domain = NA)

    for (iter in seq_len(control$iter.max)) {

        
        gradthe <- gradient(theta)     # Compute the gradient vector
        hessthe <- hessian(theta)      # Compute the Hessian matrix

        delta <- Matrix::solve(hessthe, gradthe, tol = control$tolsolve)

        if (control$trace)
            cat(iter, ': ', theta, "\n", sep = "")

        step_size <- 1L                # Initialize step size for step-halving
        

        # Step-halving loop
        while (TRUE) {
            
            new_theta <- theta - step_size * delta # Update parameter vector
            objnew_the <- objective(new_theta)

            if (control$trace)
                cat("Old, new, old - new objective:", 
                    objthe, objnew_the, objthe - objnew_the, "\n")

            # Objective function failed to be reduced or is infinite
            if (!is.finite(objnew_the) || (objnew_the > objthe + 1e-6)) {
                step_size <- step_size / 2         # Reduce the step size

                if (control$trace) 
                    cat("Step size reduced to", step_size, "\n")

                if (step_size <= control$minstepsize) {
                    msg <- paste("Step size ", step_size, 
                                 " has reduced below minstepsize")
                    return(list(par = theta, objective = objthe, convergence = 1, 
                                message = msg)) 
                }
            } else {
                theta  <- new_theta	# accept the new parameter vector
                oldobj <- objthe
                objthe <- objnew_the
                break
            }
            
        }

        
        # Convergence check - must meet 3 criteria
        if ((objthe <= oldobj + 1e-6 && (oldobj - objthe < control$objtol)) &&
            (max(abs(gradthe)) < control$gradtol) &&
            (max(abs(delta)) < control$paramtol))

            return(list(par            = theta,
                        objective      = objthe,
                        convergence    = 0,
                        message        = "Normal convergence"))
         
    }

    msg <- paste("Reached", control$iter.max, "iterations without convergence")
    return(list(par = theta, objective = objthe, convergence = 1, message = msg)) 
}



###################################################
### code chunk number 3: glm
###################################################
library("free1way.docreg")
(x <- matrix(c(10, 5, 7, 11, 8, 9), nrow = 2))
d <- expand.grid(y = relevel(gl(2, 1), "2"), t = gl(3, 1))
d$x <- c(x)
m <- glm(y ~ t, data = d, weights = x, family = binomial())
(cf <- coef(m))


###################################################
### code chunk number 4: glm-op
###################################################
F <- plogis
f <- dlogis
op <- optim(par = c("mt2" = 0, "mt3" = 0, "(Intercept)" = 0), 
            fn = .nll, gr = .nsc, 
            x = x, method = "BFGS", hessian = TRUE)
cbind(glm = c(cf[-1] * -1, cf[1]), free1way = op$par)
logLik(m)
-op$value


###################################################
### code chunk number 5: glm-H
###################################################
fp <- function(x) 
{
    p <- plogis(x)
    p * (1 - p)^2 - p^2 * (1 - p)
}
H <- .hes(op$par, x)
### analytical covariance of parameters
solve(H$Z - crossprod(H$X, Matrix::solve(H$A, H$X)))
### numerical covariance
solve(op$hessian)[1:2,1:2]
### from glm
vcov(m)[-1,-1]


###################################################
### code chunk number 6: glm-free1way
###################################################
obj <- .free1wayML(as.table(x), link = logit())
obj$coefficients
-obj$value
### analytical covariance
obj$vcov


###################################################
### code chunk number 7: glm-stratum
###################################################
(x <- as.table(array(c(10, 5, 7, 11, 8, 9,
                        9, 4, 8, 15, 5, 4), dim = c(2, 3, 2))))
d <- expand.grid(y = relevel(gl(2, 1), "2"), t = gl(3, 1), s = gl(2, 1))
d$x <- c(x)
m <- glm(y ~ 0 + s + t, data = d, weights = x, family = binomial())
logLik(m)
(cf <- coef(m))


###################################################
### code chunk number 8: glm-op-stratum
###################################################
xl <- .table2list(x)$xlist
op <- optim(par = c("mt2" = 0, "mt3" = 0, "(Intercept 1)" = 0, "(Intercept 2)" = 0), 
            fn = .snll, gr = .snsc, 
            x = xl, 
            method = "BFGS", 
            hessian = TRUE)
cbind(glm = c(cf[-(1:2)] * -1, cf[1:2]), free1way = op$par)
logLik(m)
-op$value


###################################################
### code chunk number 9: glm-H-stratum
###################################################
### analytical covariance of parameters
solve(.shes(op$par, xl))
### numerical covariance
solve(op$hessian)[1:2,1:2]
### from glm
vcov(m)[-(1:2),-(1:2)]


###################################################
### code chunk number 10: glm-free1way-strata
###################################################
obj <- .free1wayML(as.table(x), link = logit())
obj$coefficients
-obj$value
### analytical covariance
obj$vcov


###################################################
### code chunk number 11: Newton
###################################################



###################################################
### code chunk number 12: Newton-test
###################################################
N <- 10000
P <- 30
X <- matrix(rnorm(N * P), ncol = P)
y <- X %*% runif(P) + rnorm(nrow(X))

f <- function(par) sum((y - X %*% par)^2)
g <- function(par) colSums(- 2 * c(y - X %*% par) * X)
h <- function(par) 2 * crossprod(X)

start <- runif(P)

cf <- .NewtonRaphson(start = start, objective = f, gradient = g, hessian = h)

cf2 <- coef(m <- lm(y ~ 0 + X))
all.equal(sum((y - fitted(m))^2), cf$objective)
all.equal(unname(cf$par), unname(cf2))


###################################################
### code chunk number 13: workhorse
###################################################
N <- 10
a <- matrix(c(5, 6, 4,
                    3, 5, 7,
                    3, 4, 5,
                    3, 5, 6,
                    0, 0, 0,
                    4, 6, 5), ncol = 3, byrow = TRUE)
x <- as.table(array(c(a[1:3,], a[-(1:3),]), dim = c(3, 3, 2)))
x
ret <- .free1wayML(x, logit())
ret[c("value", "par")]
cf <- ret$par
cf[1:2] <- cf[1:2] + .5
### new2old parameterisation
c(cf[1:2], cf[3], log(cf[4] - cf[3]), cf[5])
### profile for cf[1:2]
.free1wayML(x, logit(), start = cf, fix = 1:2)[c("value", "par")]
### profile for cf[2]
.free1wayML(x, logit(), start = cf, fix = 2)[c("value", "par")]
### evaluate log-likelihood at cf
.free1wayML(x, logit(), start = cf, 
            fix = seq_along(ret$par))[c("value", "par")]


###################################################
### code chunk number 14: SW
###################################################
set.seed(29)
w <- gl(2, 15)
(s <- .SW(r <- rank(u <- runif(length(w))), model.matrix(~ 0 + w)))
ps <- .resample(r, model.matrix(~ 0 + w), B = 100000)
ps$testStat^2
mean(abs(ps$permStat) > abs(ps$testStat) - .Machine$double.eps)
pchisq(ps$testStat^ifelse(ps$DF == 1, 2, 1), df = ps$DF, lower.tail = FALSE)
### exactly the same
kruskal.test(u ~ w)
library("coin")
### almost the same
kruskal_test(u ~ w, distribution = approximate(100000))


###################################################
### code chunk number 15: Wexact
###################################################
wilcox_test(u ~ w, distribution = "exact")
free1way(u ~ w, exact = TRUE)


###################################################
### code chunk number 16: Wexact-le
###################################################
wilcox_test(u ~ w, distribution = "exact", alternative = "less")
print(free1way(u ~ w, exact = TRUE), alternative = "greater")


###################################################
### code chunk number 17: Wexact-gr
###################################################
wilcox_test(u ~ w, distribution = "exact", alternative = "greater")
print(free1way(u ~ w, exact = TRUE), alternative = "less")


###################################################
### code chunk number 18: formula (eval = FALSE)
###################################################
## y ~ groups | blocks


###################################################
### code chunk number 19: free
###################################################
x
### asymptotic permutation test
(ft <- free1way(x))
coef(ft)
vcov(ft)
### Wald per parameter
summary(ft)
library("multcomp")
summary(glht(ft), test = univariate())

### global Wald
summary(ft, test = "Wald")
summary(glht(ft), test = Chisqtest())

### Rao score, Permutation score, LRT
summary(ft, test = "Rao")
summary(ft, test = "Permutation")
summary(ft, test = "LRT")

### Wald confidence intervals, unadjusted
confint(glht(ft), calpha = univariate_calpha())
confint(ft, test = "Wald")

### Rao and LRT intervals
confint(ft, test = "Rao")
confint(ft, test = "LRT")


###################################################
### code chunk number 20: wilcox-formula
###################################################
N <- 25
w <- gl(2, N)
y <- rlogis(length(w), location = c(0, 1)[w])

#### link = logit is default
ft <- free1way(y ~ w)

### Wald 
summary(ft)

### Permutation test
wilcox.test(y ~ w, alternative = "greater", correct = FALSE, exact = FALSE)$p.value
pvalue(wilcox_test(y ~ w, alternative = "greater"))
summary(ft, test = "Permutation", alternative = "less")$p.value
wilcox.test(y ~ w, alternative = "less", correct = FALSE, exact = FALSE)$p.value
pvalue(wilcox_test(y ~ w, alternative = "less"))
summary(ft, test = "Permutation", alternative = "greater")$p.value
wilcox.test(y ~ w, correct = FALSE, exact = FALSE)$p.value
kruskal.test(y ~ w)$p.value
pvalue(wilcox_test(y ~ w))
summary(ft, test = "Permutation")$p.value

### Wald tests
summary(ft, test = "Wald", alternative = "less")
summary(ft, test = "Wald", alternative = "greater")
summary(ft, test = "Wald")

### Rao score tests
summary(ft, test = "Rao", alternative = "less")
summary(ft, test = "Rao", alternative = "greater")
summary(ft, test = "Rao")

### LRT (only two-sided)
summary(ft, test = "LRT")

### confidence intervals for log-odds ratios
confint(ft, test = "Permutation")
confint(ft, test = "LRT")
confint(ft, test = "Wald")
confint(ft, test = "Rao")

### confidence interval for "Wilcoxon Parameter" = PI = AUC
confint(ft, test = "Rao", what = "AUC")

### comparison with rms::orm
library("rms")
rev(coef(or <- orm(y ~ w)))[1]
coef(ft)
logLik(or)
logLik(ft)
vcov(or)[2,2]
vcov(ft)
ci <- confint(or)
ci[nrow(ci),]
confint(ft, test = "Wald")


###################################################
### code chunk number 21: mh
###################################################
example(mantelhaen.test, echo = FALSE)
mantelhaen.test(UCBAdmissions, correct = FALSE)
ft <- free1way(UCBAdmissions)
summary(ft, test = "Wald")
exp(coef(ft))
exp(confint(ft, test = "Wald"))
exp(sapply(dimnames(UCBAdmissions)[[3L]], function(dept)
       confint(free1way(UCBAdmissions[,,dept]), test = "Permutation")))
sapply(dimnames(UCBAdmissions)[[3L]], function(dept)
       fisher.test(UCBAdmissions[,,dept], conf.int = TRUE)$conf.int)


###################################################
### code chunk number 22: pt
###################################################
prop.test(UCBAdmissions[,,1], correct = FALSE)
summary(free1way(UCBAdmissions[,,1]), test = "Rao")


###################################################
### code chunk number 23: kw
###################################################
example(kruskal.test, echo = FALSE)
kruskal.test(x ~ g)
free1way(x ~ g)


###################################################
### code chunk number 24: sw
###################################################
library("survival")
N <- 10
nd <- expand.grid(g = gl(3, N), s = gl(3, N))
nd$tm <- rexp(nrow(nd))
nd$ev <- TRUE
cm <- coxph(Surv(tm, ev) ~ g + strata(s), data = nd)

(ft <- free1way(tm ~ g | s, data = nd, link = "cloglog"))
coef(cm)
coef(ft)
vcov(cm)
vcov(ft)
### Rao score tests
summary(cm)$sctest
summary(ft, test = "Rao")
### likelihood ratio tests
summary(cm)$logtest
summary(ft, test = "LRT")
### Wald tests
summary(cm)$waldtest
summary(ft, test = "Wald")
### asymptotic permutation tests
survdiff(Surv(tm, ev) ~ g + strata(s), data = nd, rho = 0)[c("chisq", "pvalue")]
summary(ft, test = "Permutation")
library("coin")
independence_test(Surv(tm, ev) ~ g | s, data = nd, ytrafo = function(...)
                  trafo(..., numeric_trafo = logrank_trafo, block = nd$s), 
                  teststat = "quad")


###################################################
### code chunk number 25: Peto
###################################################
survdiff(Surv(tm, ev) ~ g + strata(s), data = nd, rho = 1)[c("chisq", "pvalue")]
(ft <- free1way(tm ~ g | s, data = nd, link = "logit"))
summary(ft)
summary(ft, test = "Rao")
summary(ft, test = "LRT")
summary(ft, test = "Wald")
summary(ft, test = "Permutation")


###################################################
### code chunk number 26: sw
###################################################
library("survival")
data("GBSG2", package = "TH.data")
cm <- coxph(Surv(time, cens) ~ horTh + strata(tgrade), data = GBSG2)

ft <- with(GBSG2, free1way(Surv(time, cens) ~ horTh | tgrade, 
                           link = "cloglog"))
coef(cm)
coef(ft)
vcov(cm)
vcov(ft)
### Rao score tests
summary(cm)$sctest
summary(ft, test = "Rao")
### likelihood ratio tests
summary(cm)$logtest
summary(ft, test = "LRT")
### Wald tests
summary(cm)$waldtest
summary(ft, test = "Wald")
### asymptotic permutation tests
survdiff(Surv(time, cens) ~ horTh + strata(tgrade), data = GBSG2, rho = 0)[c("chisq", "pvalue")]
summary(ft, test = "Permutation")
independence_test(Surv(time, cens) ~ horTh | tgrade, data = GBSG2, ytrafo = function(...)
                  trafo(..., numeric_trafo = logrank_trafo, block = GBSG2$tgrade), 
                  teststat = "quad")


###################################################
### code chunk number 27: Peto
###################################################
survdiff(Surv(time, cens) ~ horTh + strata(tgrade), data = GBSG2, rho = 1)[c("chisq", "pvalue")]
(ft <- with(GBSG2, free1way(Surv(time, cens) ~ horTh | tgrade, 
                            link = "logit")))
summary(ft, test = "Rao")
summary(ft, test = "LRT")
summary(ft, test = "Wald")
summary(ft, test = "Permutation")


###################################################
### code chunk number 28: sw
###################################################
library("survival")
GBSG2$str <- cut(GBSG2$tsize, breaks = c(0, 1:9 * 10, Inf))
cm <- coxph(Surv(time, cens) ~ horTh + strata(str), data = GBSG2)

ft <- with(GBSG2, free1way(Surv(time, cens) ~ horTh | str, 
                           link = "cloglog"))
coef(cm)
coef(ft)
vcov(cm)
vcov(ft)
### Rao score tests
summary(cm)$sctest
summary(ft, test = "Rao")
### likelihood ratio tests
summary(cm)$logtest
summary(ft, test = "LRT")
### Wald tests
summary(cm)$waldtest
summary(ft, test = "Wald")
### asymptotic permutation tests
survdiff(Surv(time, cens) ~ horTh + strata(str), data = GBSG2, rho = 0)[c("chisq", "pvalue")]
summary(ft, test = "Permutation")


###################################################
### code chunk number 29: Peto
###################################################
survdiff(Surv(time, cens) ~ horTh + strata(str), data = GBSG2, rho = 1)[c("chisq", "pvalue")]
(ft <- with(GBSG2, free1way(Surv(time, cens) ~ horTh | str, 
                            link = "logit")))
summary(ft, test = "Rao")
summary(ft, test = "LRT")
summary(ft, test = "Wald")
summary(ft, test = "Permutation")


###################################################
### code chunk number 30: normal
###################################################
nd$y <- rnorm(nrow(nd))
free1way(y ~ g | s, data = nd, link = "probit")
independence_test(y ~ g | s, data = nd, ytrafo = function(...)
                  trafo(..., numeric_trafo = normal_trafo, block = nd$s), 
                  teststat = "quad")


###################################################
### code chunk number 31: friedman-data
###################################################
example(friedman.test, echo = FALSE)

### Myles Hollander & Wolfe (2014, Example 7.1, page 294)
boxplot(RoundingTimes, xlab = "Method", ylab = "Rounding-First-Base Time", 
        las = 1)
matplot(t(RoundingTimes), add = TRUE, type = "l", 
        lty = 1, lwd = 2, col = rgb(.1, .1, .1, .1))

me <- colnames(RoundingTimes)
d <- expand.grid(me = factor(me, labels = me, levels = me),
                 id = factor(seq_len(nrow(RoundingTimes))))
d$time <- c(t(RoundingTimes))


###################################################
### code chunk number 32: friedman
###################################################
friedman.test(RoundingTimes)
(ft <- free1way(time ~ me | id, data = d))


###################################################
### code chunk number 33: friedman-add
###################################################
summary(ft)
library("multcomp")
glht(ft, linfct = mcp(me = "Tukey"))


###################################################
### code chunk number 34: friedman-link
###################################################
logLik(ft)
logLik(free1way(time ~ me | id, data = d, link = "probit"))
logLik(free1way(time ~ me | id, data = d, link = "cloglog"))
logLik(free1way(time ~ me | id, data = d, link = "loglog"))


###################################################
### code chunk number 35: McNemar
###################################################
example(mcnemar.test, echo = FALSE)
# set-up data frame with survey outcomes for voters
s <- gl(2, 1, labels = dimnames(Performance)[[1L]])
survey <- gl(2, 1, labels = c("1st", "2nd"))
nvoters <- c(Performance)
x <- expand.grid(survey = survey, voter = factor(seq_len(sum(nvoters))))
x$performance <- c(rep(s[c(1, 1)], nvoters[1]), rep(s[c(2, 1)], nvoters[2]),
                   rep(s[c(1, 2)], nvoters[3]), rep(s[c(2, 2)], nvoters[4]))
# note that only those voters changing their minds are relevant
mcn <- free1way(xtabs(~ performance + survey + voter, data = x))
# same result as mcnemar.test w/o continuity correction
print(mcn)
# X^2 statistic
summary(mcn, test = "Permutation")$statistic^2
mcnemar.test(Performance, correct = FALSE)
# Wald inference
summary(mcn)
confint(mcn, test = "Wald")
### because the model is saturated, the link function doesn't affect the
### p-value (but the coefficients are of course different)
free1way(xtabs(~ performance + survey + voter, data = x), link = "probit")


###################################################
### code chunk number 36: incomplete
###################################################
data("taste", package = "daewr")
### highly discrete
table(taste$score)
summary(free1way(score ~ recipe | panelist, data = taste))


###################################################
### code chunk number 37: Tukey
###################################################
tk <- free1way(Ozone ~ Month, data = airquality)
library("multcomp")
confint(glht(tk, linfct = mcp(Month = "Tukey")))


###################################################
### code chunk number 38: plot-ex
###################################################
tk <- free1way(Ozone ~ Month, data = airquality)
plot(tk, las = 1)


###################################################
### code chunk number 39: plot-ex-cdf
###################################################
plot(tk, cdf = TRUE, model = FALSE, las = 1)


###################################################
### code chunk number 40: plot-ex-ecdf
###################################################
aq <- subset(airquality, !is.na(Ozone))
aq$r <- match(aq$Ozone, sort(unique(aq$Ozone)))
library("latticeExtra")
plot(ecdfplot(~ r, data = aq, groups = Month, col = 1:5))


###################################################
### code chunk number 41: ppplot
###################################################
y <- rlogis(50)
x <- rlogis(50, location = 3)
ppplot(y, x, conf.level = .95, las = 1)


###################################################
### code chunk number 42: ppplot-savage
###################################################
ppplot(y, x, conf.args = list(link = "cloglog", type = "Wald", 
                              col = NA, border = NULL),
       conf.level = .95, las = 1)


###################################################
### code chunk number 43: rfree1way
###################################################
(logOR <- c(log(1.5), log(2)))
nd <- rfree1way(150, delta = logOR)
coef(ft <- free1way(y ~ groups, data = nd))
sqrt(diag(vcov(ft)))
logLik(ft)
nd$y <- qchisq(nd$y, df = 3)
coef(ft <- free1way(y ~ groups, data = nd))
sqrt(diag(vcov(ft)))
logLik(ft)
N <- 25
pvals <- replicate(Nsim, 
{
  nd <- rfree1way(n = N, blocks = 2, delta = c(.25, .5), alloc_ratio = 2)
  summary(free1way(y ~ groups | blocks, data = nd), test = "Permutation")$p.value
})

power.free1way.test(n = N, blocks = 2, delta = c(.25, .5), alloc_ratio = 2)
mean(pvals < .05)


###################################################
### code chunk number 44: rfree1waysurv
###################################################
N <- 1000
nd <- rfree1way(N, delta = 1, link = "cloglog")
nd$C <- rfree1way(n = N, delta = 1, offset = -c(qlogis(.25), qlogis(.5)), 
                  link = "cloglog")$y
nd$y <- Surv(pmin(nd$y, nd$C), nd$y < nd$C)
### check censoring probability
1 - tapply(nd$y[,2], nd$groups, mean)
summary(free1way(y ~ groups, data = nd, link = "cloglog"))
summary(coxph(y ~ groups, data = nd))


###################################################
### code chunk number 45: power.prop.test
###################################################
delta <- log(1.5)
power.prop.test(n = 25, p1 = .5, p2 = plogis(qlogis(.5) - delta))
power.free1way.test(n = 25, prob = c(.5, .5), delta = delta)


###################################################
### code chunk number 46: power.odds.test
###################################################
prb <- matrix(c(.25, .25, .25, .25,
                .10, .20, .30, .40), ncol = 2)
colnames(prb) <- c("s1", "s2")
power.free1way.test(n = 20, prob = prb, 
                    strata_ratio = 2,
                    alloc_ratio = c(1.5, 2, 2), 
                    delta = log(c("low" = 1.25, "med" = 1.5, "high" = 1.75)))


###################################################
### code chunk number 47: wilcox
###################################################

delta <- log(3)
N <- 15
w <- gl(2, N)
pw <- numeric(Nsim)
for (i in seq_along(pw)) {
    y <- rlogis(length(w), location = c(0, delta)[w])
    pw[i] <- wilcox.test(y ~ w)$p.value
}
mean(pw < .05)

power.free1way.test(n = N, delta = delta)

### approximate formula in Hmisc::popower
library("Hmisc")
popower(p = rep(1 / N, N), odds.ratio = exp(delta), n = 2 * N)


###################################################
### code chunk number 48: kruskal
###################################################
delta <- c("B" = log(2), "C" = log(3))
N <- 15
w <- gl(3, N)
pw <- numeric(Nsim)
for (i in seq_along(pw)) {
    y <- rlogis(length(w), location = c(0, delta)[w])
    pw[i] <- kruskal.test(y ~ w)$p.value
}
mean(pw < .05)

power.free1way.test(n = N, delta = delta)


###################################################
### code chunk number 49: table
###################################################
prb <- rep.int(1, 4) / 4
pw <- numeric(Nsim)
cf <- matrix(0, nrow = Nsim, ncol = length(delta))
colnames(cf) <- names(delta)
for (i in seq_along(pw)) {
    nd <- rfree1way(n = N, prob = prb, delta = delta)
    ft <- free1way(y ~ groups, data = nd)
    cf[i,] <- coef(ft)
    pw[i] <- summary(ft, test = "Permutation")$p.value
}
mean(pw < .05)
boxplot(cf, las = 1, ylab = expression(hat(delta)))
points(c(1:2), delta, pch = 19, col = "red")
power.free1way.test(n = N, prob = prb, delta = delta)


###################################################
### code chunk number 50: stable
###################################################
prb <- cbind(S1 = rep(1, 4), 
             S2 = c(1, 2, 1, 2), 
             S3 = 1:4)
dimnames(prb) <- list(Ctrl = paste0("i", seq_len(nrow(prb))),
                      Strata = colnames(prb))

pw <- numeric(Nsim)
cf <- matrix(0, nrow = Nsim, ncol = length(delta))
colnames(cf) <- names(delta)
for (i in seq_along(pw)) {
    nd <- rfree1way(n = N, prob = prb, delta = delta)
    ft <- free1way(y ~ groups | blocks, data = nd)
    cf[i,] <- coef(ft)
    pw[i] <- summary(ft, test = "Permutation")$p.value
}
mean(pw < .05)
boxplot(cf, las = 1, ylab = expression(hat(delta)))
points(c(1:2), delta, pch = 19, col = "red")


###################################################
### code chunk number 51: powertest
###################################################
power.free1way.test(n = N, prob = prb, delta = delta, seed = 3)
power.free1way.test(power = .8, prob = prb, delta = delta, seed = 3)
power.free1way.test(n = 19, prob = prb, delta = delta, seed = 3)


###################################################
### code chunk number 52: Jeffreys
###################################################
N <- 20
w <- gl(2, N)
y <- rnorm(length(w), mean = c(-2, 3)[w])

x <- free1way(y ~ w, link = "probit")
coef(x)
logLik(x)

pll <- function(cf) {

    start <- x$par
    start[1] <- cf
    x$profile(start, fix = 1)
}

### https://doi.org/10.1111/j.0006-341X.2001.00114.x
### https://doi.org/10.1111/j.1467-9876.2012.01057.x
### https://doi.org/10.1186/s12874-017-0313-9
### https://files.osf.io/v1/resources/fet4d_v3/providers/osfstorage/682fb176db88f967facacb5a?format=pdf&action=download&direct&version=1
### https://doi.org/10.1002/sim.6537
### https://doi.org/10.1007/s11222-023-10217-3
### https://arxiv.org/abs/2510.06465
fun <- function(cf) {
    ret <- pll(cf)
    ret$value - .5 * determinant(ret$hessian, logarithm = TRUE)$modulus
}

ci <- confint(x, level = .99, test = "Wald")
grd <- seq(from = ci[1], to = ci[2], length.out = 50)

optim(coef(x), fn = fun, method = "Brent", 
      lower = min(grd), upper = max(grd))[c("par", "value")]


###################################################
### code chunk number 53: MPL_Jeffreys
###################################################
free1way(y ~ w, link = "probit", MPL_Jeffreys = TRUE)
