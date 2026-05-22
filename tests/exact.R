
pkgs <- c("free1way.docreg", "coin")
if (!all(sapply(pkgs, require, character.only = TRUE)))
    quit()

set.seed(29)

N <- 15
w <- gl(2, N)[-(1:floor(N/2))]

pvals <- clwr <- cupr <- matrix(0, nrow = 10, ncol = 2)

for (i in seq_len(nrow(pvals))) {
    y <- round(rlogis(length(w), location = c(0, .5)[w]))
    fe <- free1way(y ~ w, exact = TRUE)
    fa <- free1way(y ~ w, B = 100000)
    pvals[i,1] <- summary(fe, test = "Perm")$p.value
    pvals[i,2] <- summary(fa, test = "Perm")$p.value
    ci <- confint(fe)
    clwr[i,1] <- ci[1]
    cupr[i,1] <- ci[2]
    ci <- confint(fa)
    clwr[i,2] <- ci[1]
    cupr[i,2] <- ci[2]
}

max(abs(pvals[,1] - pvals[,2])) < .005
max(abs(clwr[,1] - clwr[,2])) < .05
max(abs(cupr[,1] - cupr[,2])) < .05

