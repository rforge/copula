## Copyright (C) 2010 Marius Hofert and Martin Maechler
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

stopifnot(require(Runuran), require(optimx))

## ==== GIG generator and related functions ====================================

## generator (parameters theta_1 in IR, theta_2 in (0,Inf))
## note: theta_1!=0 and theta_2=0 reduces to copClayton@psi(t, 1/theta_1)
psi.GIG <- function(t, theta)
    (1+t)^(-theta[1]/2)*besselK(theta[2]*sqrt(1+t),nu=theta[1])/
    besselK(theta[2],nu=theta[1])

## generator inverse 
## note: the initial interval is only valid for theta_1 > -1/2
##       it can be large, but it's finite; it is based on an inequality about 
##       modified Bessel functions of the third kind (besselK) given in 
##       Paris (1984) ("An inequality for the Bessel function J_v(vx)", 
##       SIAM J. Math. Anal. 15, 203--205)
psiInv.GIG <- function(t, theta, upper=if(theta[1] > -0.5) function(x) (1-log(x)/theta[2])^2-1 else 
                       stop("initial interval for psiInv.GIG fails"), ...){ 
    res <- numeric(length(t))
    res[t0 <- t==0] <- Inf
    res[t1 <- t==1] <- 0
    t. <- t[!(t0 | t1)]
    l <- length(t.)
    up <- upper(t.)
    res <- unlist(lapply(1:l, function(i){
	uniroot(function(t..) psi.GIG(t.., theta)-t.[i],
                interval=c(0, up[i]), ...)$root
    }))
    if(is.matrix(t)) matrix(res, ncol=ncol(t)) else res
}

## generator derivatives
psiDabs.GIG <- function(t, theta, degree=1, n.MC=0, log=FALSE){
    if(n.MC>0){
        V <- V0.GIG(n.MC, theta)
        lx <- -V %*% t(t) + degree*log(V)
        res <- colMeans(exp(lx))
        r <- if(log) log(res) else res
    }else{
        res <- numeric(n <- length(t))
        res[is0 <- t == 0] <- Inf
        res[isInf <- is.infinite(t)] <- -Inf
        if(any(n0Inf <- !(is0 | isInf))){
            res[n0Inf] <- degree*log(theta[2]/2)-((theta[1]+degree)/2)*log1p(t) +
                log(besselK(theta[2]*sqrt(1+t), nu=theta[1]+degree, expon.scaled=TRUE)) - 
                    log(besselK(theta[2], nu=theta[1], expon.scaled=TRUE)) -
                        (sqrt(1+t)-1)*theta[2]
        }
        r <- if(log) res else exp(res)
    }
    if(is.matrix(t)) matrix(r, ncol=ncol(t)) else r
}

## derivative of the generator inverse
psiInvD1abs.GIG <- function(t, theta, log=FALSE){
    res <- -psiDabs.GIG(psiInv.GIG(t, theta), theta, log=TRUE)
    if(log) res else exp(res)
}

## derivative of the generator inverse with given psiInv = psiInv.GIG(t, theta)
psiInvD1absInv.GIG <- function(psiInv, theta, log=FALSE){
    res <- -psiDabs.GIG(psiInv, theta, log=TRUE)
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
    psiI <- psiInv.GIG(u., theta=theta) # this is costly if d is large
    psiI.sum <- rowSums(psiI)
    ## main part
    if(n.MC > 0){ # Monte Carlo
        res[n01] <- psiDabs.GIG(psiI.sum, theta=theta, degree=d, n.MC=n.MC, log=TRUE) +
            rowSums(psiInvD1absInv.GIG(psiI, theta=theta, log=TRUE))
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
    if(log) log(2*ud(dens, 2*x)) else 2*ud(dens, 2*x)
}

## Kendall's tau for fixed theta_1 as a function in theta_2 [vectorized]
## FIXME: this is numerically difficult, take, e.g., tau.GIG(c(1,400))
##        or even tau.GIG(c(0.00001,0.00001))
tau.GIG <- function(theta, ...){
    if(!is.matrix(theta)) theta <- rbind(theta)
    integrand <- function(t, theta)
        t/(1+t)^(theta[1]+1)*(besselK(theta[2]*sqrt(1+t), nu=theta[1]+1)/
                              besselK(theta[2], nu=theta[1]))^2
    ## t/(1+t)^(theta[1]+1)*exp(2*(log(besselK(theta[2]*sqrt(1+t), nu=theta[1]+1))-log(besselK(theta[2], nu=theta[1]))))
    res <- 1-theta[,2]^2*apply(theta, 1, function(x) 
                               integrate(integrand, lower=0, upper=Inf, theta=x, ...)$value)
    names(res) <- NULL
    res
}

## tau inverse for all parameters fixed except one (given by NA)
## note: initial interval has to be specified [e.g., c(1e-12,100)] since there 
##       are no adequate bounds known
tauInv.GIG <- function(tau, theta, interval, ...){
    stopifnot(length(i <- which(is.na(theta))) == 1)
    replace(theta, i, uniroot(function(th) tau.GIG(replace(theta, i, th))-tau, 
                              interval=interval, ...)$root)
}

## generate vectors of random variates from a GIG copula
rnacopula.GIG <- function(n, d, theta)
    psi.GIG(matrix(rexp(n*d), nrow=n, ncol=d)/V0.GIG(n, theta), theta)

