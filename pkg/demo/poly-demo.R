library(nacopula)
library(animation)
library(lattice)

options(warn=1)

## ==== plot of poly* (= polyG, polyJ) for all methods =========================

## ==== expected evaluation points for estimating Gumbel and Joe copulas ====

eep.fun <- function(family, alpha, d, n.MC=5000){
    vapply(alpha, function(alph)
       {
           th <- 1/alph
           cop <- onacopulaL(family, list(th, 1:d))
           U <- rnacopula(n.MC, cop)
           switch(family,
                  Gumbel =
              {
                  mean(rowSums(cop@copula@psiInv(U, th))^alph)
              },
                  Joe =
              {
                  U. <- (1-U)^th
                  lh <- rowSums(log1p(-U.)) # log(h(..))
                  l1_h <- log(-expm1(lh))
                  mean(exp(lh-l1_h))
              }, stop("wrong family in eep.fun()"))
       }, NA_real_)
}

## ==== plot function for a vector of alphas ====

plot.poly <- function(family, xlim, ylim, method, alpha, d){
    len <- length(alpha)
    cols <- colorRampPalette(c("red", "orange", "darkgreen", "turquoise", "blue"),
                             space="Lab")(len)
    switch(family,
           "Gumbel" = {
               fun <- nacopula:::polyG
               str <- "G"
           },
           "Joe" = {
               fun <- nacopula:::polyJ
               str <- "J"
           },
       {stop("wrong family in plot.poly")})
    for(k in 1:len){
        curve(fun(log(x), alpha=alpha[k], d=d, method=method, log=TRUE),
              from=xlim[1], to=xlim[2], main=paste("poly",str,
                                        "(log(x), alpha=..., d=",d,
                                        ", method=",method,", log=TRUE)",sep=""),
              xlab="x", ylab=paste("log(poly",str,"(log(x), ...))", sep=""),
              add=(k>1), lwd=1.4, col=cols[k], ylim=ylim)
    }
    label <- as.expression(lapply(1:len, function(i)
        substitute(alpha==a, list(a=alpha[i]))))
    legend("bottomright", label, bty="n", lwd=1.4, col=cols)
}

## ==== animation in alpha ====

## animation of poly* functions
## m = number of frames
## d = dimension
## method = method for polyG
poly.ani <- function(family, m, d, method, xlim, ylim){
    switch(family,
           "Gumbel" = {
               fun <- nacopula:::polyG
               str <- "G"
           },
           "Joe" = {
               fun <- nacopula:::polyJ
               str <- "J"
           },
       {stop("wrong family in plot.poly")})
    alphas <- (1:m)/(m+1) # alphas
    eep <- eep.fun(family, alphas, d) # corresponding expected evaluation points for Gumbel
    x <- seq(xlim[1], xlim[2], length.out=1000)
    lx <- log(x)
    res <- lapply(1:m, function(i){
        if(i %% 5 == 1) print(paste(round(i/m*100),"% done",sep="")) # progress
        y <- fun(lx, alpha=alphas[i], d=d, method=method,
                 log=TRUE)
        p <- xyplot(y~x, type="l", xlab="x",
                    ylab=paste("log(poly",str,"(log(x), ...))",sep=""), aspect=1,
                    xlim=xlim, ylim=ylim, key=list(x=0.35, y=0.1, lines=list(lty=1,
                                                                  col="black"),
                                          text=list(paste("expected x-value for alpha=",
                                          alphas[i],sep=""))),
                    panel=function(...){
                        panel.xyplot(...)
                        panel.abline(v=eep[i]) # vertical line at expected value
                    }, main=paste("poly",str,"(log(x), alpha=",alphas[i],
                       ", d=",d,", method=",method,", log=TRUE)",sep=""))
        list(y=y, plot=p)
    })
    res
}

## ==== Gumbel =================================================================

## ==== plots for small d ====

family <- "Gumbel"
alpha <- c(0.99, 0.5, 0.01) # alphas; plot largest first, so that all values are visible
d <- 5
xlim <- c(1e-16, 1000)
ylim <- c(-40, 40)
(ev <- eep.fun(family, alpha=alpha, d=d)) # 4.927077 2.462882 1.025498
stopifnot(all(xlim[1] < ev, ev < xlim[2]))

meths <- eval(formals(nacopula:::polyG)$method)
for(i in seq_along(meths)){
    plot.poly(family, xlim=xlim, ylim=ylim, method=meths[i], alpha=alpha, d=d)
    Sys.sleep(2)
}
## => all are fine, even for the much larger range than the expected values

## ==== plots for large d ====

d <- 100
xlim <- c(1e-16, 200)
ylim <- c(300, 600)
(ev <- eep.fun(family, alpha, d)) # 96.135490 11.128448  1.060105
stopifnot(all(xlim[1] < ev, ev < xlim[2]))

## method == "pois"
plot.poly(family, xlim=xlim, ylim=ylim, method="pois", alpha=alpha, d=d)
## => problems for small and moderate alpha

## method == "pois.direct"
plot.poly(family, xlim=xlim, ylim=ylim, method="pois.direct", alpha=alpha, d=d)
## => problems for small and moderate alpha

## method == "binomial.coeff"
plot.poly(family, xlim=xlim, ylim=ylim, method="binomial.coeff", alpha=alpha, d=d)
## => same as "pois"

## method == "stirling"
plot.poly(family, xlim=xlim, ylim=ylim, method="stirling", alpha=alpha, d=d)
## => problems only for large alphas

## method == "stirling.horner"
plot.poly(family, xlim=xlim, ylim=ylim, method="stirling.horner", alpha=alpha, d=d)
## => same as "stirling"

## other methods
plot.poly(family, xlim=xlim, ylim=ylim, method="sort", alpha=alpha, d=d) # log(< 0)
plot.poly(family, xlim=xlim, ylim=ylim, method="horner", alpha=alpha, d=d) # log(< 0)
plot.poly(family, xlim=xlim, ylim=ylim, method="direct", alpha=alpha, d=d) # log(< 0)
plot.poly(family, xlim=xlim, ylim=ylim, method="dJoe", alpha=alpha, d=d) # okay for large alpha

## ==== run time comparison of the methods that worked for some parameter ====

set.seed(1)
x <- runif(100000, min=0.01, max=120)
lx <- log(x)

## pois: for large alpha (where it works)
system.time(y.pois <- nacopula:::polyG(lx, alpha=0.99, d=d, method="pois",
                                       log=TRUE))[[1]]
## => 8.91s
stopifnot(all(!is.nan(y.pois))) # check

## pois.direct: for large alpha (where it works)
system.time(y.pois.direct <- nacopula:::polyG(lx, alpha=0.99, d=d,
                                              method="pois.direct", log=TRUE))[[1]]
## => 6.80s
stopifnot(all(!is.nan(y.pois.direct))) # check

## binomial.coeff: for large alpha (where it works)
system.time(y.binomial.coeff <- nacopula:::polyG(lx, alpha=0.99, d=d,
                                                 method="binomial.coeff",
                                                 log=TRUE))[[1]]
## => 8.72s
stopifnot(all(!is.nan(y.binomial.coeff))) # check

## stirling: for moderate alpha (where it works)
system.time(y.stirling <- nacopula:::polyG(lx, alpha=0.5, d=d,
                                           method="stirling", log=TRUE))[[1]]
## => 1.92s
stopifnot(all(!is.nan(y.stirling))) # check

## stirling.horner: for moderate alpha (where it works)
system.time(y.stirling.horner <- nacopula:::polyG(lx, alpha=0.5, d=d,
                                                  method="stirling.horner",
                                                  log=TRUE))[[1]]
## => 2.79s
stopifnot(all(!is.nan(y.stirling.horner))) # check

## dJoe: for large alpha (where it works)
system.time(y.dJoe <- nacopula:::polyG(lx, alpha=0.99, d=d, method="dJoe",
                                       log=TRUE))[[1]]
## => 2.28s
stopifnot(all(!is.nan(y.dJoe))) # check

## conclusion:
## - fastest for large alpha: "dJoe", "pois.direct"
## - fastest for small and moderate alpha: "stirling"
## - further methods tried: pulling out max() for "stirling" => does not increase precision

## ==== more detailed graphical precision comparison in d = 100 ====

## dJoe
m <- 49
ylim <- c(200, 700)
polyG.ani.dJoe <- poly.ani(family, m, d=d, method="dJoe", xlim=c(1e-16,200),
                           ylim=ylim)
saveHTML(for(i in 1:m) print(polyG.ani.dJoe[[i]]$plot))
## => works for alpha >= 0.75

## pois.direct
polyG.ani.pois.direct <- poly.ani(family, m, d=d, method="pois.direct",
                                  xlim=c(1e-16,200), ylim=ylim)
saveHTML(for(i in 1:m) print(polyG.ani.pois.direct[[i]]$plot))
## => works for the whole range of *expected* values, esp. for alpha >= 0.72

## stirling
polyG.ani.stirling <- poly.ani(family, m, d=d, method="stirling",
                               xlim=c(1e-16,200), ylim=ylim)
saveHTML(for(i in 1:m) print(polyG.ani.stirling[[i]]$plot))
## => works for alpha <= 0.56

## ==== check default method ====

## comparison with Maple (Digits = 100)
v1 <- nacopula:::polyG(log(1), alpha=0.01, d=d, log=TRUE)
v2 <- nacopula:::polyG(log(1), alpha=0.5, d=d, log=TRUE)
v3 <- nacopula:::polyG(log(1), alpha=0.99, d=d, log=TRUE)
stopifnot(all.equal(c(v1,v2,v3), c(354.52779560, 356.56733266, 350.99662083))) # comparison with Maple

v1 <- nacopula:::polyG(log(17), alpha=0.01, d=d, log=TRUE)
v2 <- nacopula:::polyG(log(17), alpha=0.5, d=d, log=TRUE)
v3 <- nacopula:::polyG(log(17), alpha=0.99, d=d, log=TRUE)
stopifnot(all.equal(c(v1,v2,v3), c(358.15179523, 374.67231305, 370.20372192))) # comparison with Maple

v1 <- nacopula:::polyG(log(77), alpha=0.01, d=d, log=TRUE)
v2 <- nacopula:::polyG(log(77), alpha=0.5, d=d, log=TRUE)
v3 <- nacopula:::polyG(log(77), alpha=0.99, d=d, log=TRUE)
stopifnot(all.equal(c(v1,v2,v3), c(362.38428102, 422.83827969, 435.36899283), tol=1e-6)) # comparison with Maple

## animation in alpha
polyG.ani.default <- poly.ani(family, m, d=d, method="default",
                              xlim=c(1e-16,200), ylim=ylim)
saveHTML(for(i in 1:m) print(polyG.ani.default[[i]]$plot))


## ==== Joe ====================================================================

## ==== plots for small d ====

set.seed(1)
family <- "Joe"
alpha <- c(0.05, 0.5, 0.99) # alphas; plot smallest first, so that all values are visible
d <- 5
xlim <- c(1e-16, 1e120)
ylim <- c(0, 1200)
(ev <- eep.fun(family, alpha, d, n.MC=100000)) # 5.618664e+79 6.923153e+03 5.912850e-02; varies a lot for different runs!
if(!all(xlim[1] < ev, ev < xlim[2])) warning("ev outside xlim")

meths <- eval(formals(nacopula:::polyJ)$method)
for(i in seq_along(meths)){
    plot.poly(family, xlim=xlim, ylim=ylim, method=meths[i], alpha=alpha, d=d)
    Sys.sleep(2)
}
## => does not work for any reasonable x range

## ==== plots for large d ====

set.seed(1)
d <- 100
xlim <- c(1e-16, 1e120)
ylim <- c(0, 30000)
(ev <- eep.fun(family, alpha, d, n.MC=100000)) # 2.119430e+78 8.011466e+02 2.040599e-04; varies a lot for different runs!
if(!all(xlim[1] < ev, ev < xlim[2])) warning("ev outside xlim")

## method == "log.poly"
plot.poly(family, xlim=xlim, ylim=ylim, method="log.poly", alpha=alpha, d=d)
## => flawlessly (as far as visible)

## method == "log1p"
plot.poly(family, xlim=xlim, ylim=ylim, method="log1p", alpha=alpha, d=d)
## => flawlessly (as far as visible)

## method == "poly"
plot.poly(family, xlim=xlim, ylim=ylim, method="poly", alpha=alpha, d=d)
## => does not work for any reasonable x range

## ==== run time comparison of the methods that worked for some parameter ====

set.seed(1)
x <- runif(100000, min=0.01, max=1e100)
lx <- log(x)

## log.poly:
system.time(y.log.poly <- nacopula:::polyJ(lx, alpha=0.5, d=d, method="log.poly",
                                           log=TRUE))[[1]]
## => 2.701s
stopifnot(all(!is.nan(y.log.poly))) # check

## log1p:
system.time(y.log1p <- nacopula:::polyJ(lx, alpha=0.5, d=d,
                                        method="log1p", log=TRUE))[[1]]
## => 4.118s
stopifnot(all(!is.nan(y.log1p))) # check

## conclusion: use log.poly as default

## comparison with Maple (Digits = 100)
v1 <- nacopula:::polyJ(log(1), alpha=0.01, d=d, log=TRUE)
v2 <- nacopula:::polyJ(log(1), alpha=0.5, d=d, log=TRUE)
v3 <- nacopula:::polyJ(log(1), alpha=0.99, d=d, log=TRUE)
stopifnot(all.equal(c(v1,v2,v3), c(395.73694325, 393.08027226, 386.96715831))) # comparison with Maple

v1 <- nacopula:::polyJ(log(1e20), alpha=0.01, d=d, log=TRUE)
v2 <- nacopula:::polyJ(log(1e20), alpha=0.5, d=d, log=TRUE)
v3 <- nacopula:::polyJ(log(1e20), alpha=0.99, d=d, log=TRUE)
stopifnot(all.equal(c(v1,v2,v3), c(4918.2008336, 4915.3815020, 4909.1039909))) # comparison with Maple

v1 <- nacopula:::polyJ(log(1e100), alpha=0.01, d=d, log=TRUE)
v2 <- nacopula:::polyJ(log(1e100), alpha=0.5, d=d, log=TRUE)
v3 <- nacopula:::polyJ(log(1e100), alpha=0.99, d=d, log=TRUE)
stopifnot(all.equal(c(v1,v2,v3), c(23154.67477009, 23151.85543852, 23145.57792740))) # comparison with Maple

## animation in alpha
polyJ.ani.default <- poly.ani(family, m, d=d, method="log.poly",
                              xlim=c(1e-16,1e120), ylim=ylim)
saveHTML(for(i in 1:m) print(polyJ.ani.default[[i]]$plot))
## => rather extreme but seems to be fine

