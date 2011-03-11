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
psiInv.GIG <- function(t, theta, upper=NULL, ...){ 
    if(is.null(upper)){
	if(theta[1] > -0.5) function(x) (1-log(x)/theta[2])^2-1 else 
        stop("initial interval for psiInv.GIG fails")
    }
    res <- numeric(length(t))
    is0 <- t==0
    is1 <- t==1
    n01 <- !(is0 | is1)
    res[is0] <- Inf
    res[is1] <- 0
    t. <- t[n01]
    up <- upper(t.)
    res[n01] <- unlist(lapply(1:length(t.), function(i){
	uniroot(function(t..) psi.GIG(t.., theta)-t.[i],
                interval=c(0, up[i]), ...)$root
    }))
    if(is.matrix(t)) matrix(res, ncol=ncol(t)) else res
}

## generator derivatives
psiDabs.GIG <- function(t, theta, degree=1, n.MC=0, log=FALSE){
    res <- numeric(length(t))
    iInf <- is.infinite(t)
    res[iInf] <- -Inf # log(0)
    if(any(!iInf)){
	t. <- t[!iInf]
        if(n.MC>0){
            V <- V0.GIG(n.MC, theta)
            res[!iInf] <- nacopula:::lsum(-V %*% t(t.) + degree*log(V) - log(n.MC))
        }else{
            res[!iInf] <- degree*log(theta[2]/2)-((theta[1]+degree)/2)*log1p(t.) +
                log(besselK(theta[2]*sqrt(1+t.), nu=theta[1]+degree, expon.scaled=TRUE)) - 
                    log(besselK(theta[2], nu=theta[1], expon.scaled=TRUE)) -
                        (sqrt(1+t.)-1)*theta[2]
        }
    }
    r <- if(log) res else exp(res)
    if(is.matrix(t)) matrix(r, ncol=ncol(t)) else r
}

## derivative of the generator inverse
psiInvD1abs.GIG <- function(t, theta, log=FALSE){
    res <- -psiDabs.GIG(psiInv.GIG(t, theta=theta), theta=theta, log=TRUE)
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
        res[n01] <- psiDabs.GIG(psiI.sum, theta=theta, degree=d, n.MC=n.MC, log=TRUE) -
            rowSums(psiDabs.GIG(psiI, theta=theta, log=TRUE))
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
## note: theta can be a vector of the form (nu,...,nu,theta,...,theta)
## theta.min denots the smallest theta for which no limiting formula is used
tau.GIG <- function(theta, theta.min=1e-4, ...){
    if(!is.matrix(theta)) theta <- matrix(theta, ncol=2)
    lim <- theta[,1]>0 & theta[,2]<theta.min # indices for which limiting formula for theta=0 is used (since numerically challenging)
    res <- numeric(nrow(theta))
    res[lim] <- 1/(1+2*theta[lim,1])
    if(any(!lim)){
	integrand <- function(t, theta){
            st <- sqrt(1+t)
            t*(theta[2]*besselK(theta[2]*st, nu=theta[1]+1)/(st^(theta[1]+1)*
                                             besselK(theta[2], nu=theta[1])))^2
        }
	res[!lim] <- 1-apply(theta[!lim,,drop=FALSE], 1, function(x) 
                             integrate(integrand, lower=0, upper=Inf, theta=x, ...)$value)
    }
    names(res) <- NULL
    res
}

## inverse of Kendall's tau for all parameters fixed except the one given as NA
## note: initial interval has to be specified [e.g., c(1e-30,100)] since there 
##       are no adequate bounds known
tauInv.GIG <- function(tau, theta, interval, iargs=list(), theta.min=1e-4, ...){
    stopifnot(length(i <- which(is.na(theta))) == 1)
    if(i!=1){ # theta[1]=nu is given
	if(theta[1]<=0) stop("tauInv.GIG: only implemented for positive nu") 
	max.tau <- 1/(1+2*theta[1])
	if(tau>max.tau) stop(paste("tauInv.GIG: maximum attainable tau for nu=",
                                   theta[1]," is ",round(max.tau,8),sep="")) # this assumes tau to be falling in theta (numerical experiments show this behavior)
    }else{ # theta[2]=theta given
	if(theta[2]<theta.min) return(2*tau/(1-tau)) # limiting formula
    }
    replace(theta, i, uniroot(function(th) do.call(tau.GIG, c(list(replace(theta, 
                                                                           i, th), 
                                                                   theta.min=theta.min), 
                                                              iargs))-tau, 
                              interval=interval, ...)$root)
}

## generate vectors of random variates from a GIG copula
rnacopula.GIG <- function(n, d, theta)
    psi.GIG(matrix(rexp(n*d), nrow=n, ncol=d)/V0.GIG(n, theta), theta)

