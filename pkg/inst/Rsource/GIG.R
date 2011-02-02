stopifnot(require(Runuran), require(optimx))

## ==== GIG generator and related functions ====================================

## generator (parameters theta_1 in IR, theta_2 in (0,Inf))
## note: theta_1!=0 and theta_2=0 reduces to copClayton@psi(t, 1/theta_1)
psi.GIG <- function(t, theta)
    (1+t)^(-theta[1]/2)*besselK(theta[2]*sqrt(1+t),theta[1])/besselK(theta[2],theta[1])

## generator inverse 
## note: the initial interval is only valid for theta_1 > -1/2
##       it can be large, but it's finite; it is based on an inequality about 
##       modified Bessel functions of the third kind (besselK) given in 
##       Paris (1984) ("An inequality for the Bessel function J_v(vx)", 
##       SIAM J. Math. Anal. 15, 203--205)
psiInv.GIG <- function(t, theta, ...){ 
    if(theta[1] <= -0.5) stop("initial interval for psiInv.GIG fails")
    res <- numeric(length(t))
    res[t==0] <- Inf
    res[t==1] <- 0
    unlist(lapply(t, function(t.){
	uniroot(function(t..) psi.GIG(t..,theta)-t.,
                interval=c(0, (1-log(t.)/theta[2])^2-1), ...)$root
    }))
}

## generator derivatives
psiDabs.GIG <- function(t, theta, degree=1, n.MC=0, log=FALSE){
    if(n.MC>0){
        V0. <- V0.GIG(n.MC, theta)
        l.V0. <- degree*log(V0.)
        lx <- -V0. %*% t(t) + l.V0.
	res <- colMeans(exp(lx))
        if(log) log(res) else res
    }else{
        res <- numeric(n <- length(t))
        res[is0 <- t == 0] <- Inf
        res[isInf <- is.infinite(t)] <- -Inf
        if(any(n0Inf <- !(is0 | isInf))){
            res[n0Inf] <- d*log(theta[2]/2)-(sqrt(1+t)-1)*theta[2]+
                log(besselK(theta[2]*sqrt(1+t),theta[1],expon.scaled=TRUE))-
                    log(besselK(theta[2],theta[1],expon.scaled=TRUE))-((theta[1]+d)/2)*log1p(t)
        }
        if(log) res else exp(res)
    }
}

## derivative of the generator inverse
psiInvD1abs.GIG <- function(t, theta, log=FALSE){
    res <- -psiDabs.GIG(psiInv.GIG(t, theta), theta, log=TRUE)
    if(log) res else exp(res)
}

## density of the GIG copula
dacopula.GIG <- function(u, theta, n.MC=0, log=FALSE){
    if(!is.matrix(u)) u <- rbind(u)
    if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
    ## f() := NaN outside and on the boundary of the unit hypercube
    res <- rep.int(NaN, n <- nrow(u))
    n01 <- apply(u,1,function(x) all(0 < x, x < 1)) # indices for which density has to be evaluated
    if(!any(n01)) return(res)
    ## auxiliary results
    u. <- u[n01,, drop=FALSE]
    psiI <- matrix(psiInv.GIG(u., theta=theta), ncol=d) # this is the costly part if d is large
    psiI.sum <- rowSums(psiI)
    ## main part
    if(n.MC > 0){ # Monte Carlo
	res[n01] <- psiDabs.GIG(psiI.sum, theta=theta, degree=d, n.MC=n.MC, log=TRUE) -
            rowSums(matrix(psiInvD1abs.GIG(u., theta=theta, log=TRUE), ncol=d))
    }else{ # explicit
        ## auxiliary function for density evaluation
        h <- function(t, k, theta, log=FALSE){
            s1pt <- sqrt(1+t)
            B <- besselK(theta[2]*s1pt, nu=theta[1]+k)
            if(log) log(B)-(theta[1]+k)*log(s1pt) else B/s1pt^(theta[1]+k)
        }
        ## result
        res[n01] <- h(psiI.sum, k=d, theta=theta, log=TRUE) + (d-1)*
            log(besselK(theta[2], nu=theta[1])) - rowSums(h(psiI, k=1, theta=theta, 
                                  log=TRUE))
    }
    if(log) res else exp(res)
}

## V0 random number generator
V0.GIG <- function(n, theta){
    dens <- udgig(theta[1], 1, theta[2]^2) # the three args lambda, psi, chi for the GIG density as on page 497 in McNeil, Frey, Embrechts (2005)
    gen <- pinvd.new(dens) # works fast, via approximation of the quantile function by piecewise polynomials
    ur(gen, n)/2
}

## density of V0
dV0.GIG <- function(x, theta, log=FALSE){
    dens <- udgig(theta[1], 1, theta[2]^2)
    if(log) log(ud(dens, x)) else ud(dens, x)
}

## Kendall's tau for fixed theta_1 as a function in theta_2
## note: this is numerically difficult, take, e.g., tau.GIG(1, 400)
##       or even tau.GIG(0.00001,0.00001)
tau.GIG <- function(theta, method=c("direct", "transformation"), ...){
    method <- match.arg(method)
    switch(method,
           "direct" = {
               integrand <- function(t, theta){
                   exp(log(t)-(theta[1]+1)*log1p(t)-2*sqrt(1+t)*theta[2]+
                       2*log(besselK(theta[2]*sqrt(1+t), theta[1]+1, expon.scaled=TRUE)))
               }
               1-(theta[2]/besselK(theta[2], theta[1]))^2 *
                   integrate(integrand, lower=0, upper=Inf, theta=theta, ...)$value	
           },
           "transformation" = { # note: numerically not more stable...
               integrand <- function(t, theta)
                   (t^2-theta[2]^2)/t^(4*theta[1]+3)*(t^(theta[1]+1)*besselK(t, theta[1]+1, expon.scaled=TRUE)*exp(-t))^2
               1-2*theta[2]^(4*theta[1])*(theta[2]^theta[1]*besselK(theta[2], theta[1], expon.scaled=TRUE)*exp(-theta[2]))^(-2)*
                   integrate(integrand, lower=theta[2], upper=Inf, theta=theta, ...)$value
           },
       {stop(sprintf("unsupported method '%s' in psiDabsMC", method))})  
}

## vectorized wrapper of tau.GIG in theta_2 for given theta_1 [can be a number or a vector]
tau.GIG.vec <- function(theta1, theta2, method=c("direct", "transformation"), ...){
    n.1 <- length(theta1)
    n.2 <- length(theta2)
    stopifnot(n.1==1 || n.1==n.2)
    if(n.1==1) theta1 <- rep(theta1, n.2)
    theta <- cbind(theta1, theta2)
    apply(theta, 1, tau.GIG, method=method, ...)
}

## tau inverse for fixed theta_1
## note: initial interval has to be specified [e.g., c(1e-12,100)] since there 
##       are no adequate bounds known
tauInv.GIG <- function(tau, theta1, interval, ...)
    uniroot(function(th) tau.GIG(c(theta1,th))-tau, interval=interval, ...)$root

## generate vectors of random variates from a GIG copula
rnacopula.GIG <- function(n, d, theta)
    psi.GIG(matrix(rexp(n*d), nrow=n, ncol=d)/V0.GIG(n, theta), theta)

