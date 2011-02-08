library(nacopula)
options(warn=1)

## ==== plot of polyG for all methods ==========================================

## ==== plot function for a vector of alphas ====

plot.polyG <- function(from, to, method, alpha, d){
    len <- length(alpha)
    cols <- colorRampPalette(c("red", "orange", "darkgreen", "turquoise", "blue"), 
                             space="Lab")(len)
    for(k in 1:len){ 
        curve(nacopula:::polyG(log(x), alpha=alpha[k], d=d, method=method, log=TRUE),
              from=from, to=to, main=paste("polyG(log(x), alpha=..., d=",d,
                                ", method=",method,", log=TRUE)",sep=""), 
              xlab="x", ylab=paste("log(polyG(log(x), ...))", sep=""), 
              add=(k>1), lwd=1.4, col=cols[k])
    }	
    label <- as.expression(lapply(1:len, function(i) 
        substitute(alpha==a, list(a=alpha[i]))))
    legend("bottomright", label, bty="n", lwd=1.4, col=cols)
}

## ==== expected evaluation points for estimating a Gumbel copula ==== 

eep.fun <- function(alpha, d, n.MC=5000){
    len <- length(alpha)
    ex <- numeric(len) # expected x
    for(k in 1:len){
	th <- 1/alpha[k]
	cop <- onacopulaL("Gumbel", list(th, 1:d))
	U <- rnacopula(n.MC, cop)
	ex[k] <- mean(rowSums(copGumbel@psiInv(U, th))^alpha[k])
    }
    ex
}

## ==== plots for small d ====

alpha <- c(0.99, 0.5, 0.01) # alphas; plot largest first, so that all values are visible
d <- 5
from <- 1e-16
to <- 1000
(ev <- eep.fun(alpha,d)) # 4.927077 2.462882 1.025498
stopifnot(all(from < ev, ev < to))
meths <- eval(formals(nacopula:::polyG)$method)
for(i in seq_along(meths)){
    plot.polyG(from, to, method=meths[i], alpha=alpha, d=d)
    Sys.sleep(2.4)
}
## => all are fine, even for the much larger range than the expected values

## ==== plots for large d ====

alpha <- c(0.99, 0.5, 0.01) # alphas; plot largest first, so that all values are visible
d <- 100
from <- 1e-16
to <- 200
(ev <- eep.fun(alpha,d)) # 95.911722 11.247274 1.053746
stopifnot(all(from < ev, ev < to))

## method == "pois"
plot.polyG(from, to, method="pois", alpha=alpha, d=d)
## => problems for small and moderate alpha 

## method == "pois.direct"
plot.polyG(from, to, method="pois.direct", alpha=alpha, d=d)
## => problems for small and moderate alpha

## method == "binomial.coeff"
plot.polyG(from, to, method="binomial.coeff", alpha=alpha, d=d)
## => same as "pois"

## method == "stirling"
plot.polyG(from, to, method="stirling", alpha=alpha, d=d)
## => problems only for large alphas

## method == "stirling.horner"
plot.polyG(from, to, method="stirling.horner", alpha=alpha, d=d)
## => same as "stirling"

## other methods
plot.polyG(from, to, method="sort", alpha=alpha, d=d) # log(< 0)
plot.polyG(from, to, method="horner", alpha=alpha, d=d) # log(< 0)
plot.polyG(from, to, method="direct", alpha=alpha, d=d) # log(< 0)
plot.polyG(from, to, method="dJoe", alpha=alpha, d=d) # okay for large alpha

## ==== run time comparison of the methods that worked ====

set.seed(1)
x <- runif(100000, min=0.01, max=120) 
lx <- log(x)

## pois: for large alpha (where it works)
system.time(y.pois <- nacopula:::polyG(lx, alpha=0.99, d=100, method="pois", 
                                       log=TRUE))[[1]]
## => 8.91s
stopifnot(all(!is.nan(y.pois))) # check

## pois.direct: for large alpha (where it works)
system.time(y.pois.direct <- nacopula:::polyG(lx, alpha=0.99, d=100, 
                                              method="pois.direct", log=TRUE))[[1]]
## => 6.80s
stopifnot(all(!is.nan(y.pois.direct))) # check

## binomial.coeff: for large alpha (where it works)
system.time(y.binomial.coeff <- nacopula:::polyG(lx, alpha=0.99, d=100, 
                                                 method="binomial.coeff", 
                                                 log=TRUE))[[1]]
## => 8.72s
stopifnot(all(!is.nan(y.binomial.coeff))) # check

## stirling: for moderate alpha (where it works)
system.time(y.stirling <- nacopula:::polyG(lx, alpha=0.5, d=100, 
                                           method="stirling", log=TRUE))[[1]]
## => 1.92s
stopifnot(all(!is.nan(y.stirling))) # check

## stirling.horner: for moderate alpha (where it works)
system.time(y.stirling.horner <- nacopula:::polyG(lx, alpha=0.5, d=100, 
                                                  method="stirling.horner", 
                                                  log=TRUE))[[1]]
## => 2.79s
stopifnot(all(!is.nan(y.stirling.horner))) # check

## dJoe: for large alpha (where it works)
system.time(y.dJoe <- nacopula:::polyG(lx, alpha=0.99, d=100, method="dJoe", 
                                       log=TRUE))[[1]]
## => 2.28s
stopifnot(all(!is.nan(y.dJoe))) # check

## conclusion: 
## fastest for large alpha: "dJoe", "pois.direct"
## fastest for small and moderate alpha: "stirling"
## further methods tried: pulling out max() for "stirling" => does not increase precision

## ==== more detailed graphical precision comparison in d = 100 ====

library(animation)
library(lattice)

## animation of polyG
## m = number of frames
## d = dimension
## method = method for polyG
polyG.animate <- function(m, d, method, xlim=c(1e-16,200), ylim){
    alphas <- (1:m)/(m+1) # alphas
    eep <- eep.fun(alphas, d) # corresponding expected evaluation points for Gumbel
    x <- seq(xlim[1],xlim[2],length.out=1000)
    lx <- log(x)
    res <- lapply(1:m, function(i){
        if(i %% 5 == 1) print(paste(round(i/m*100),"% done",sep="")) # progress
        y <- nacopula:::polyG(lx, alpha=alphas[i], d=d, method=method, 
                              log=TRUE)
        p <- xyplot(y~x, type="l", xlab="x", 
                    ylab="log(polyG(log(x), ...))", aspect=1, 
                    xlim=xlim, ylim=ylim, key=list(x=0.35, y=0.1, lines=list(lty=1, 
                                                                  col="black"), 
                                          text=list(paste("expected x-value for alpha=",
                                          alphas[i],sep=""))),
                    panel=function(...){
                        panel.xyplot(...)
                        panel.abline(v=eep[i]) # vertical line at expected value
                    }, main=paste("polyG(log(x), alpha=",alphas[i],
                       ", d=",d,", method=",method,", log=TRUE)",sep=""))
        list(y=y, plot=p)
    })
    res
}

## dJoe
m <- 49
polyG.animate.dJoe <- polyG.animate(m, d=100, method="dJoe", ylim=c(200,700))
saveHTML(for(i in 1:m) print(polyG.animate.dJoe[[i]]$plot))
## => works for alpha >= 0.75

## pois.direct
polyG.animate.pois.direct <- polyG.animate(m, d=100, method="pois.direct", 
                                           ylim=c(200,700))
saveHTML(for(i in 1:m) print(polyG.animate.pois.direct[[i]]$plot))
## => works for the whole range of *expected* values, esp. for alpha >= 0.72

## stirling
polyG.animate.stirling <- polyG.animate(m, d=100, method="stirling", 
                                        ylim=c(200,700))
saveHTML(for(i in 1:m) print(polyG.animate.stirling[[i]]$plot))
## => works for alpha <= 0.56

## ==== combine the methods ====

## comparison with Maple (Digits = 100)
nacopula:::polyG(log(1), alpha=0.01, d=100, log=TRUE)  # Maple: 354.52779560
nacopula:::polyG(log(1), alpha=0.5, d=100, log=TRUE)   # Maple: 356.56733266
nacopula:::polyG(log(1), alpha=0.99, d=100, log=TRUE)  # Maple: 350.99662083

nacopula:::polyG(log(17), alpha=0.01, d=100, log=TRUE) # Maple: 358.15179523
nacopula:::polyG(log(17), alpha=0.5, d=100, log=TRUE)  # Maple: 374.67231305
nacopula:::polyG(log(17), alpha=0.99, d=100, log=TRUE) # Maple: 370.20372192

nacopula:::polyG(log(77), alpha=0.01, d=100, log=TRUE) # Maple: 362.38428102
nacopula:::polyG(log(77), alpha=0.5, d=100, log=TRUE)  # Maple: 422.83827969
nacopula:::polyG(log(77), alpha=0.99, d=100, log=TRUE) # Maple: 435.36899283
# => great

## animate in alpha
polyG.animate.default <- polyG.animate(m, d=100, method="default", ylim=c(200,700))
saveHTML(for(i in 1:m) print(polyG.animate.default[[i]]$plot))

