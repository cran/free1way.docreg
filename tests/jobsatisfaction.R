
pkgs <- c("free1way.docreg", "coin", "multcomp")
if (!all(sapply(pkgs, require, character.only = TRUE)))
    quit()

jstab <- jobsatisfaction

js <- as.data.frame(jstab)
js$Job.Satisfaction <- ordered(js$Job.Satisfaction,
      levels = dimnames(jstab)[[2]],
      labels = dimnames(jstab)[[2]])
js$jsnum <- unclass(js$Job.Satisfaction)

jsall <- js[rep(1:nrow(js), js$Freq),]

### all the same
free1way(Job.Satisfaction ~ Income | Gender, data = jsall)
free1way(jsnum ~ Income | Gender, data = jsall)
kruskal_test(jsnum ~ Income | Gender, data = jsall)
free1way(Job.Satisfaction ~ Income | Gender, data = js, weights = Freq)
free1way(jsnum ~ Income | Gender, data = js, weights = Freq)
free1way(aperm(jstab, perm = c(2, 1, 3)))

### w/o coin
example(mantelhaen.test, echo = FALSE)
(ft <- free1way(aperm(Satisfaction, perm = c(2, 1, 3))))

### post-hoc
glht(ft, linfct = mcp(Income = "Dunnett"))
