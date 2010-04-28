library(nacopula)

set.seed(101)
X <- rstable1(1e4, alpha=1, beta=1, gamma= .25, delta=1)
summary(X)
   ##   Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
   ## 0.4235    0.8959    1.1400    2.7670    1.6140 2302.0000 

dX <- density(log(X))
tit <- expression( density(log( rstable1(n, alpha == 1, beta == 1, gamma == .25, delta == 1))) )
plot(dX, main= tit, xlim = c(-1.5, 5), ylim = c(0, 1.1))
for(i in 1:20)
    lines(density(log(X2 <- rstable1(1e4, alpha=1, beta=1, gamma= .25, delta=1))), col=2)
