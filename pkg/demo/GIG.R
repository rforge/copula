library(nacopula) # for (nested) Archimedean copulas and auxiliary functions
library(Runuran) # for generating random numbers of a GIG distribution etc.
library(optimx) # for optimization 
options(warn=1) # print warnings as they occur

#### ==== Two-parameter Archimedean GIG copula =================================

## ==== GIG generator and related functions ====================================

## generator (parameters nu in IR, theta in (0,Inf))
## note: nu!=0 and theta=0 reduces to copClayton@psi(t, 1/nu)
psi.GIG <- function(t, nu, theta)
    (1+t)^(-nu/2)*besselK(theta*sqrt(1+t),nu)/besselK(theta,nu)

## generator inverse 
## note: the initial interval is only valid for nu > -1/2
##       it can be large, but it's finite; it is based on an inequality about 
##       modified Bessel functions of the third kind (besselK) given in 
##       Paris (1984) ("An inequality for the Bessel function J_v(vx)", 
##       SIAM J. Math. Anal. 15, 203--205)
psiInv.GIG <- function(t, nu, theta, ...){ 
    if(nu <= -0.5) stop("initial interval for psiInv.GIG fails")
    res <- numeric(length(t))
    res[t==0] <- Inf
    res[t==1] <- 0
    unlist(lapply(t, function(t.){
	uniroot(function(t..) psi.GIG(t..,nu,theta)-t.,
                interval=c(0, (1-log(t.)/theta)^2-1), ...)$root
    }))
}

## generator derivatives
psiDabs.GIG <- function(t, nu, theta, degree=1, n.MC=0, log=FALSE){
    if(n.MC>0){
        V0. <- V0.GIG(n.MC, nu, theta)
        l.V0. <- degree*log(V0.)
        lx <- -V0. %*% t(t) + l.V0.
	res <- colMeans(exp(lx))
        if(log) log(res) else res
    }else{
        res <- numeric(n <- length(t))
        res[is0 <- t == 0] <- Inf
        res[isInf <- is.infinite(t)] <- -Inf
        if(any(n0Inf <- !(is0 | isInf))){
            res[n0Inf] <- d*log(theta/2)-(sqrt(1+t)-1)*theta+
                log(besselK(theta*sqrt(1+t),nu,expon.scaled=TRUE))-
                    log(besselK(theta,nu,expon.scaled=TRUE))-((nu+d)/2)*log1p(t)
        }
        if(log) res else exp(res)
    }
}

## derivative of the generator inverse
psiInvD1abs.GIG <- function(t, nu, theta, log=FALSE){
    res <- -psiDabs.GIG(psiInv.GIG(t, nu, theta), nu, theta, log=TRUE)
    if(log) res else exp(res)
}

## density of the GIG copula
dacopula.GIG <- function(u, nu, theta, n.MC=0, log=FALSE){
    if(!is.matrix(u)) u <- rbind(u)
    if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
    ## f() := NaN outside and on the boundary of the unit hypercube
    res <- rep.int(NaN, n <- nrow(u))
    n01 <- apply(u,1,function(x) all(0 < x, x < 1)) # indices for which density has to be evaluated
    if(!any(n01)) return(res)
    ## auxiliary results
    u. <- u[n01,, drop=FALSE]
    psiI <- matrix(psiInv.GIG(u., nu=nu, theta=theta), ncol=d) # this is the costly part if d is large
    psiI.sum <- rowSums(psiI)
    ## main part
    if(n.MC > 0){ # Monte Carlo
	res[n01] <- psiDabs.GIG(psiI.sum, nu=nu, theta=theta, degree=d, n.MC=n.MC, log=TRUE) -
            rowSums(matrix(psiInvD1abs.GIG(u., nu=nu, theta=theta, log=TRUE), ncol=d))
    }else{ # explicit
        ## auxiliary function for density evaluation
        h <- function(t, k, nu, theta, log=FALSE){
            s1pt <- sqrt(1+t)
            B <- besselK(theta*s1pt, nu=nu+k)
            if(log) log(B)-(nu+k)*log(s1pt) else B/s1pt^(nu+k)
        }
        ## result
        res[n01] <- h(psiI.sum, k=d, nu=nu, theta=theta, log=TRUE) + (d-1)*
            log(besselK(theta, nu=nu)) - rowSums(h(psiI, k=1, nu=nu, theta=theta, 
                               log=TRUE))
    }
    if(log) res else exp(res)
}

## V0 random number generator
V0.GIG <- function(n, nu, theta){
    dens <- udgig(nu, 1, theta^2) # the three args lambda, psi, chi for the GIG density as on page 497 in McNeil, Frey, Embrechts (2005)
    gen <- pinvd.new(dens) # works fast, via approximation of the quantile function by piecewise polynomials
    ur(gen, n)/2
}

## density of V0
dV0.GIG <- function(x, nu, theta, log=FALSE){
    dens <- udgig(nu, 1, theta^2)
    if(log) log(ud(dens, x)) else ud(dens, x)
}

## Kendall's tau for fixed nu as a function in theta
## note: this is numerically difficult, take, e.g., tau.GIG(1, 400)
##       or even tau.GIG(0.00001,0.00001)
tau.GIG. <- function(nu, theta, method=c("direct", "transformation"), ...){
    method <- match.arg(method)
    switch(method,
           "direct" = {
               integrand <- function(t, nu, theta){
                   exp(log(t)-(nu+1)*log1p(t)-2*sqrt(1+t)*theta+
                       2*log(besselK(theta*sqrt(1+t), nu+1, expon.scaled=TRUE)))
               }
               1-(theta/besselK(theta, nu))^2 *
                   integrate(integrand, lower=0, upper=Inf, nu=nu, theta=theta, ...)$value	
           },
           "transformation" = { # note: numerically not more stable...
               integrand <- function(t, nu, theta)
                   (t^2-theta^2)/t^(4*nu+3)*(t^(nu+1)*besselK(t, nu+1, expon.scaled=TRUE)*exp(-t))^2
               1-2*theta^(4*nu)*(theta^nu*besselK(theta, nu, expon.scaled=TRUE)*exp(-theta))^(-2)*
                   integrate(integrand, lower=theta, upper=Inf, nu=nu, theta=theta, ...)$value
           },
       {stop(sprintf("unsupported method '%s' in psiDabsMC", method))})  
}
tau.GIG <- Vectorize(tau.GIG., "theta") # vectorize in theta

## tau inverse for fixed nu
## note: initial interval has to be specified [e.g., c(1e-12,100)] since there 
##       are no adequate bounds known
tauInv.GIG <- function(tau, nu, interval, ...)
    uniroot(function(th) tau.GIG(nu, th)-tau, interval=interval, ...)$root

## lower tail dependence coefficient 
lambdaL = function(nu, theta){
    n <- length(nu)
    stopifnot(n == length(theta))
    rep(0, n)
}
lambdaLInv = function(lambda) {
    if(any(lambda != 0))
        stop("Any parameter for a GIG copula gives lambdaL = 0")
    NA * lambda
}
## upper tail dependence coefficient 
lambdaU = function(nu, theta){
    n <- length(nu)
    stopifnot(n == length(theta))
    rep(0, n)
}
lambdaUInv = function(lambda) {
    if(any(lambda != 0))
        stop("Any parameter for a GIG copula gives lambdaU = 0")
    NA * lambda
}

## ==== play-around with the GIG family ========================================

set.seed(1) # set seed

## ==== specify parameters ====

d <- 5 # dimension
tau <- 0.5 # specify Kendall's tau
nu <- 0.01 # fix the parameter nu of the GIG generator
(theta <- tauInv.GIG(tau, nu, interval=c(1e-12,100))) # compute theta such that Kendall's tau equals tau

## initial interval bounds for optimization
h.nu <- 0.05 # determines initial interval for nu
h.th <- 0.05 # determines initial interval for theta
I.nu <- c(nu-h.nu, nu+h.nu) # initial interval for nu
I.th <- c(theta-h.th, theta+h.th) # initial interval for theta

## ==== sample the copula ====

rnacopula.GIG <- function(n, d, theta, nu)
    psi.GIG(matrix(rexp(n*d), nrow=n, ncol=d)/V0.GIG(n, nu, theta), nu, theta)
n <- 1000
U <- rnacopula.GIG(n, d, theta, nu)
tau.hat <- cor(U, method="kendall")
tau.hats <- tau.hat[upper.tri(tau.hat)]
stopifnot(all.equal(rep(tau, d*(d-1)/2), tau.hats, tol=0.03)) # check tau

## plot the sample
if(FALSE) pairs(U, gap=0)

## ==== estimate the parameters ====

n <- 100 # choose small sample size due to run time
U <- rnacopula.GIG(n, d, theta, nu)

## -sum(log(density)) with x = c(nu, theta) as required for optimx
optim.dacopula.GIG <- function(x, u, n.MC=0)
    -sum(dacopula.GIG(u, x[1], x[2], n.MC=n.MC, log=TRUE))

## choose random initial point in a reasonable interval
start <- c(runif(1, min=I.nu[1], max=I.nu[2]), runif(1, min=I.th[1], max=I.th[2]))
## call optimizer for MLE
system.time(res.MLE <- optimx(par=start, fn=optim.dacopula.GIG, u=U, n.MC=FALSE, 
                              lower=c(I.nu[1],I.th[1]), upper=c(I.nu[2],I.th[2]),
                              method="bobyqa"))
res.MLE

## call optimizer for SMLE
system.time(res.SMLE <- optimx(par=start, fn=optim.dacopula.GIG, u=U, n.MC=10000, 
                               lower=c(I.nu[1],I.th[1]), upper=c(I.nu[2],I.th[2]),
                               method="bobyqa"))
res.SMLE

## note: run time is mainly determined by the evaluations of psiInv.GIG

## ==== plot Kendall's tau as a function in theta for different nu's ====

if(FALSE){
    nus <- c(0, 0.001, 0.5, 1, 5, 10)
    cols <- colorRampPalette(c("red", "orange", "darkgreen", "turquoise", "blue"), 
                             space="Lab")(length(nus))
    for(i in seq_along(nus))
        curve(tau.GIG(nu=nus[i], theta=x), 1e-12, 2, 
              main="Kendall's tau for the GIG family", ylim=c(0,1),
              xlab=expression(theta), ylab=expression(tau(nu,theta)), add=(i>1), 
              lwd=1.4, col=cols[i])
    legend("topright", paste("nu =",nus), bty="n", lwd=1.4, col=cols)
    ## conclusion: - largest range of tau's for nu = 0
    ##             - not possible to evaluate for nu < 0
    ##             - tau.GIG is numerically critical for nu > 0 close to 0             
}
