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

#### List of supported Archimedean copulas

## FIXME: Not "nice" that the names of the copula objects must be be used
##	inside some of their components.
## Possible solution:
##	use *same* environment for  (tau, tauInv, paraConstr, nestConstr)

## NOTA BENE:  Write psi(), tau(), ... functions such that they *vectorize*
## ---------   *and* do work even for (non-numeric) NULL argument
##          	{now checked in "acopula" validityMethod -- see ./AllClass.R }

### ==== Ali-Mikhail-Haq, see Nelsen (2007) p. 116, # 3 ========================

copAMH <-
    new("acopula", name = "AMH",
        ## generator
        psi = function(t,theta) { (1-theta)/(exp(t+0)-theta) },
        psiInv = function(t,theta) { log((1-theta*(1-t))/t) },
	## absolute value of generator derivatives
	psiDabs = function(t, theta, degree = 1, n.MC, log = FALSE){
	    if(!(missing(n.MC) || is.null(n.MC))){
		psiDabsMC(t,"AMH",theta,degree,n.MC,log)
            }else{
		if(theta == 0) if(log) return(-t) else return(exp(-t)) # independence
		## Note: psiDabs(0, ...) is correct
		t.et <- theta*exp(-t)
		if(log)
		    polylog(t.et, s = -degree, method = "neg", log=TRUE)+ log1p(-theta)-log(theta)
		else
		    polylog(t.et, s = -degree, method = "neg")* (1-theta)/theta
            }
        },
        ## derivatives of the generator inverse
	psiInvD1abs = function(t, theta, log = FALSE){
            if(log){
                log1p(-theta)-log(t)-log1p(-theta*(1-t))
            }else{
                (1-theta)/(t*(1-theta*(1-t)))
            }
	},
	## density
	dacopula = function(u, theta, n.MC, log = FALSE){
	    if(!is.matrix(u)) u <- rbind(u)
            if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
            ## f() := NaN outside and on the boundary of the unit hypercube
	    res <- rep.int(NaN, n <- nrow(u))
            n01 <- apply(u,1,function(x) all(0 < x, x < 1)) # indices for which density has to be evaluated
            if(any(n01)) { # if there are any u's inside the unit hypercube
                if(theta == 0){ res[n01] <- if(log) 0 else 1; return(res) } # independence
		## auxiliary results
	        u. <- u[n01,, drop=FALSE]
                u.. <- -theta*(1-u.)
                sum. <- rowSums(log(u.))
                sum.. <- rowSums(log1p(u..))
                ## main part
                if(!(missing(n.MC) || is.null(n.MC))){ # Monte Carlo
                    V <- copAMH@V0(n.MC, theta)
                    l <- d*log((1-theta)*V)
                    ln01 <- sum(n01)
		    one.u <- function(i) exp(l + (V-1)*sum.[i] - (V+1)*sum..[i])
                    res[n01] <- colMeans(matrix(unlist(lapply(1:ln01, one.u)),
                                                ncol = ln01))
                    if(log) log(res) else res
                }else{ # explicit
                    Li.arg <- theta*apply(u./(1+u..), 1, prod)
                    Li. <- polylog(Li.arg, s = -d, method = "neg", log=TRUE)
                    res[n01] <- (d+1)*log1p(-theta)-log(theta)+Li.-sum.-sum..
                    if(log) res else exp(res)
                }
	    } else res # all NA
        },
        ## parameter interval
        paraInterval = interval("[0,1)"),
        ## parameter subinterval (meant to be for *robust* optimization, root-finding etc.)
        paraSubInterval = num2interval(c(0,1-1e-12)), # 1-1e-12 corresponds to tau = 0.3333333
        ## nesting constraint
        nestConstr = function(theta0,theta1) {
            copAMH@paraConstr(theta0) &&
            copAMH@paraConstr(theta1) && theta1 >= theta0
        },
        ## V0 with density dV0 and V01 with density dV01 corresponding to
        ## LS^{-1}[exp(-V_0psi_0^{-1}(psi_1(t)))]
        V0 = function(n,theta) rgeom(n, 1-theta) + 1,
        dV0 = function(x,theta,log = FALSE) dgeom(x-1, 1-theta, log),
        V01 = function(V0,theta0,theta1) {
            rnbinom(length(V0),V0,(1-theta1)/(1-theta0))+V0
        },
        dV01 = function(x,V0,theta0,theta1,log = FALSE){
            stopifnot(length(V0) == 1 || length(x) == length(V0))
            dnbinom(x-V0,V0,(1-theta1)/(1-theta0),log = log)
        },
        ## Kendall's tau
        tau = tauAMH, ##-> ./aux-acopula.R
        ## function(th)  1 - 2*((1-th)*(1-th)*log(1-th)+th)/(3*th*th)
        ## but numerically stable, including, theta -> 0
        tauInv = function(tau, tol = .Machine$double.eps^0.25, ...) {
            if(any(tau > 1/3))
                stop("Impossible for AMH copula to attain a Kendall's tau larger than 1/3")
            sapply(tau,function(tau) {
                r <- safeUroot(function(th) tauAMH(th) - tau,
                               interval = as.numeric(copAMH@paraSubInterval),
                               Sig = +1, tol = tol, check.conv=TRUE, ...)
                r$root
            })
        },
        ## lower tail dependence coefficient lambda_l
        lambdaL = function(theta) { 0*theta },
        lambdaLInv = function(lambda) {
            if(any(lambda != 0))
                stop("Any parameter for an Ali-Mikhail-Haq copula gives lambdaL = 0")
            NA * lambda
        },
        ## upper tail dependence coefficient lambda_u
        lambdaU = function(theta) { 0*theta },
        lambdaUInv = function(lambda) {
            if(any(lambda != 0))
                stop("Any parameter for an Ali-Mikhail-Haq copula gives lambdaU = 0")
            NA * lambda
        })

### ==== Clayton, see Nelsen (2007) p. 116, #1 (slightly simpler form) =========

copClayton <-
    new("acopula", name = "Clayton",
        ## generator
        psi = function(t,theta) { (1+t)^(-1/theta) },
        psiInv = function(t,theta) { t^(-theta) - 1 },
        ## absolute value of generator derivatives
        psiDabs = function(t, theta, degree = 1, n.MC, log = FALSE){
            if(!(missing(n.MC) || is.null(n.MC))){
                psiDabsMC(t,"Clayton",theta,degree,n.MC,log)
            }else{
                ## Note: psiDabs(0, ...) is correct
                alpha <- 1/theta
                res <- lgamma(degree+alpha)-lgamma(alpha)-(degree+alpha)*log1p(t)
                if(log) res else exp(res)
            }
        },
        ## derivatives of the generator inverse
        psiInvD1abs = function(t, theta, log = FALSE){
            if(log) log(theta)-(1+theta)*log(t) else theta*t^(-(1+theta))
        },
	## density
	dacopula = function(u, theta, n.MC, log = FALSE){
	    if(!is.matrix(u)) u <- rbind(u)
	    if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
            ## f() := NaN outside and on the boundary of the unit hypercube
	    res <- rep.int(NaN, n <- nrow(u))
            n01 <- apply(u,1,function(x) all(0 < x, x < 1)) # indices for which density has to be evaluated
            if(any(n01)) { # if there are any u's inside the unit hypercube
		## auxiliary results
	        u. <- u[n01,, drop=FALSE]
                l.u <- rowSums(-log(u.))
	        ## main part
	        if(!(missing(n.MC) || is.null(n.MC))){ # Monte Carlo
	            V <- copClayton@V0(n.MC, theta)
	            l <- d*log(theta*V)
                    theta. <- 1 + theta
                    psiI.sum <- rowSums(copClayton@psiInv(u., theta))
                    ln01 <- sum(n01)
		    one.u <- function(i) exp(l + theta.*l.u[i] - V*psiI.sum[i])
	            res[n01] <- colMeans(matrix(unlist(lapply(1:ln01, one.u)),
                                                ncol = ln01))
	            if(log) log(res) else res
	        }else{ # explicit
                    alpha <- 1/theta
	            d.a <- d + alpha
                    arg <- rowSums(u.^(-theta)-1)
		    res[n01] <- lgamma(d.a)-lgamma(alpha)+ d*log(theta) +
			(1+theta)*l.u - (d.a)*log1p(arg)
	            if(log) res else exp(res)
	        }
	    } else res # all NA
	},
        ## parameter interval
        paraInterval = interval("(0,Inf)"),
        ## parameter subinterval (meant to be for *robust* optimization, root-finding etc.)
        paraSubInterval = num2interval(c(1e-12, 100)), # 100 corresponds to tau = 0.98
        ## nesting constraint
        nestConstr = function(theta0,theta1) {
            copClayton@paraConstr(theta0) &&
            copClayton@paraConstr(theta1) && theta1 >= theta0
        },
        ## V0 with density dV0 and V01 with density dV01 corresponding to
        ## LS^{-1}[exp(-V_0psi_0^{-1}(psi_1(t)))]
        V0 = function(n,theta) { rgamma(n, shape = 1/theta) },
        dV0 = function(x,theta,log = FALSE) dgamma(x, shape = 1/theta, log),
        V01 = function(V0,theta0,theta1) { retstable(alpha=theta0/theta1, V0) },
        dV01 = function(x,V0,theta0,theta1,log = FALSE){
            stopifnot(length(V0) == 1 || length(x) == length(V0))
            alpha <- theta0/theta1
            beta <- 1
            gamma <- (cos(pi/2*alpha)*V0)^(1/alpha)
            delta <- V0*(alpha == 1)
            if(log){
                V0-x+mapply(dstable,x,alpha = alpha,beta = beta,gamma = gamma,
                            delta = delta,pm = 1,log = TRUE)
            }else{
                exp(V0-x)*mapply(dstable,x,alpha = alpha,beta = beta,gamma = gamma,
                                 delta = delta,pm = 1)
            }
        },
        ## Kendall's tau
        tau = function(theta) { theta/(theta+2) },
        tauInv = function(tau) { 2*tau/(1-tau) },
        ## lower tail dependence coefficient lambda_l
        lambdaL = function(theta) { 2^(-1/theta) },
        lambdaLInv = function(lambda) { -1/log2(lambda) },
        ## upper tail dependence coefficient lambda_u
        lambdaU = function(theta) { 0*theta },
        lambdaUInv = function(lambda) {
            if(any(lambda != 0))
                stop("Any parameter for a CLayton copula gives lambdaU = 0")
            NA * lambda
        })

### ==== Frank, see Nelsen (2007) p. 116, # 5 ==================================

##' Frank object
copFrank <-
    new("acopula", name = "Frank",
        ## generator
        psi = function(t,theta) {
            -log1p(expm1(-theta)*exp(0-t))/theta
            ## == -log(1-(1-exp(-theta))*exp(-t))/theta
        },
        psiInv = function(t,theta) {
            -log(expm1(-theta*t)/expm1(-theta))
            ## == -log((exp(-theta*t)-1)/(exp(-theta)-1))
        },
        ## absolute value of generator derivatives
        psiDabs = function(t, theta, degree = 1, n.MC, log = FALSE){
            if(!(missing(n.MC) || is.null(n.MC))){
                psiDabsMC(t,"Frank",theta,degree,n.MC,log)
            }else{
                ## Note: psiDabs(0, ...) is correct
                p <- -expm1(-theta)
		if(log)
		    polylog(p*exp(-t), s = -(degree-1), method = "neg", log=TRUE) - log(theta)
		else
		    polylog(p*exp(-t), s = -(degree-1), method = "neg")/theta
            }
        },
        ## derivatives of the generator inverse
        psiInvD1abs = function(t, theta, log = FALSE){
            if(log) log(theta)-log(expm1(theta*t)) else theta/expm1(theta*t)
        },
	## density
	dacopula = function(u, theta, n.MC, log = FALSE){
	    if(!is.matrix(u)) u <- rbind(u)
	    if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
            ## f() := NaN outside and on the boundary of the unit hypercube
	    res <- rep.int(NaN, n <- nrow(u))
	    n01 <- apply(u,1,function(x) all(0 < x, x < 1)) # indices for which density has to be evaluated
	    if(any(n01)){ # if there are any u's inside the unit hypercube
		## auxiliary results
		u. <- u[n01,, drop=FALSE]
                u.sum <- rowSums(u.)
                lp <- log1p(-exp(-theta)) # log(p), p = 1-exp(-theta)
                lpu <- log1p(-exp(-theta*u.)) # log(1 - exp(-theta * u))
                lu <- rowSums(lpu)
	        ## main part
	        if(!(missing(n.MC) || is.null(n.MC))){ # Monte Carlo
	            V <- copFrank@V0(n.MC, theta)
	            l <- d*(log(theta*V)-V*lp)
	            rs <- -theta*u.sum
                    ln01 <- sum(n01)
		    one.u <- function(i) exp(rs[i] + l + (V-1)*lu[i])
	            res[n01] <- colMeans(matrix(unlist(lapply(1:ln01, one.u)),
                                                ncol = ln01))
	            if(log) log(res) else res
	        }else{ # explicit
	            Li.arg <- -expm1(-theta) * apply(exp(lpu-lp), 1, prod) ##<<-- FIXME faster
	            Li. <- polylog(Li.arg, s = -(d-1), method = "neg", log=TRUE)
	            res[n01] <- (d-1)*log(theta) + Li. - theta*u.sum - lu
	            if(log) res else exp(res)
	        }
	    } else res # all NA
	},
        ## parameter interval
        paraInterval = interval("(0,Inf)"),
        ## parameter subinterval (meant to be for *robust* optimization, root-finding etc.)
        paraSubInterval = num2interval(c(1e-12, 198)), # 198 corresponds to tau = 0.98
        ## nesting constraint
        nestConstr = function(theta0,theta1) {
            copFrank@paraConstr(theta0) &&
            copFrank@paraConstr(theta1) && theta1 >= theta0
        },
        ## V0 (algorithm of Kemp (1981)) with density dV0 and V01 with density
        ## dV01 corresponding to LS^{-1}[exp(-V_0psi_0^{-1}(psi_1(t)))]
        V0 = function(n,theta) rlog(n, -expm1(-theta), exp(-theta)),
        dV0 = function(x,theta,log = FALSE){
            if(any(x != (x <- floor(x + 0.5)))) warning("x must be integer; is rounded with floor(x+0.5) otherwise")
            if(log){
                x*log1p(-exp(-theta))-log(x*theta)
            }else{
                p <- 1-exp(-theta)
                p^x/(x*theta)
            }
        },
        V01 = function(V0, theta0, theta1, rej = 1, approx = 10000) {
            ## rej method switch: if V0*theta_0*p0^(V0-1) > rej a rejection
            ## from F01 of Joe is applied (otherwise, the sum is
            ## sampled via a logarithmic envelope for the summands)
            ## approx is the largest number of summands before asymptotics is used (see copJoe@V01)
            ## FIXME: optimal value of rej (for approx = 10000) is not clear yet; rej = 1 is not bad, however
            ##        lgammacor gets underflow warnings for rej < 1
            rF01Frank(V0, theta0, theta1, rej, approx)
        },
        dV01 = function(x,V0,theta0,theta1,log = FALSE){
            stopifnot(length(V0) == 1 || length(x) == length(V0), all(x >= V0))
            alpha <- theta0/theta1
            e0 <- exp(-theta0)
            e1 <- exp(-theta1)
            lfactor <- x*log1p(-e1)-V0*log1p(-e0)
            ljoe <- copJoe@dV01(x,V0,theta0,theta1,TRUE)
            res <- lfactor+ljoe
            if(log) res else exp(res)
        },
        ## Kendall's tau; debye_1() is from package 'gsl' :
        tau = function(theta){
            if((l <- length(theta)) == 0) return(numeric(0)) # to work with NULL
            r <- numeric(l)
            r[na <- is.na(theta)] <- NA
            r[!na] <- 1 + 4*(debye_1(theta[!na]) - 1)/theta[!na]
            r
        },
        tauInv = function(tau, tol = .Machine$double.eps^0.25, ...) {
            sapply(tau, function(tau) {
                r <- safeUroot(function(th) copFrank@tau(th) - tau,
                               interval = as.numeric(copFrank@paraSubInterval),
                               Sig = +1, tol=tol,
                               check.conv=TRUE, ...)
                r$root
            })
        },
        ## lower tail dependence coefficient lambda_l
        lambdaL = function(theta) { 0*theta },
        lambdaLInv = function(lambda) {
            if(any(lambda != 0))
                stop("Any parameter for a Frank copula gives lambdaL = 0")
            NA * lambda
        },
        ## upper tail dependence coefficient lambda_u
        lambdaU = function(theta) { 0*theta },
        lambdaUInv = function(lambda) {
            if(any(lambda != 0))
                stop("Any parameter for a Frank copula gives lambdaU = 0")
            NA * lambda
        })

### ==== Gumbel, see Nelsen (2007) p. 116, # 4 =================================

copGumbel <-
    new("acopula", name = "Gumbel",
        ## generator
        psi = function(t,theta) { exp(-t^(1/theta)) },
        psiInv = function(t,theta) { (-log(t+0))^theta },
        ## absolute value of generator derivatives
        psiDabs = function(t, theta, degree = 1, n.MC, log = FALSE){
            if(!(missing(n.MC) || is.null(n.MC))){
                psiDabsMC(t,"Gumbel",theta,degree,n.MC,log)
            }else{
                if(theta == 1) return(if(log) -t else exp(-t)) # independence
                res <- numeric(n <- length(t))
                res[is0 <- t == 0] <- Inf
                res[isInf <- is.infinite(t)] <- -Inf
                if(any(n0Inf <- !(is0 | isInf))) {
                    t. <- t[n0Inf]
                    t.a <- t. ^ (alpha <- 1/theta)
                    sum. <- psiDabsGpoly(t.a, alpha, degree)
                    res[n0Inf] <- -t.a - degree*log(t.) + log(sum.)
                }
                if(log) res else exp(res)
            }
        },
        ## derivatives of the generator inverse
        psiInvD1abs = function(t, theta, log = FALSE){
            if(log){
                l.t <- log(t)
                log(theta)+(theta-1)*log(-l.t)-l.t
            }else{
                theta*(-log(t))^(theta-1)/t
            }
        },
	## density
	dacopula = function(u, theta, n.MC, log = FALSE){
	    if(!is.matrix(u)) u <- rbind(u)
	    if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
            ## f() := NaN outside and on the boundary of the unit hypercube
	    res <- rep.int(NaN, n <- nrow(u))
	    n01 <- apply(u,1,function(x) all(0 < x, x < 1)) # indices for which density has to be evaluated
	    if(any(n01)){ # if there are any u's inside the unit hypercube
	        if(theta == 1){ res[n01] <- if(log) 0 else 1; return(res) } # independence
		## auxiliary results
	        u. <- u[n01,, drop=FALSE]
                psiI. <- rowSums(copGumbel@psiInv(u., theta))
                l.u <- -log(u.)
                ll.u <- log(l.u)
	        ## main part
	        if(!(missing(n.MC) || is.null(n.MC))){ # Monte Carlo
	            V <- copGumbel@V0(n.MC, theta)
	            l <- d*log(theta*V)
	            sum. <- rowSums((theta-1)*ll.u + l.u)
                    ln01 <- sum(n01)
		    one.u <- function(i) exp(l - V*psiI.[i] + sum.[i])
	            res[n01] <- colMeans(matrix(unlist(lapply(1:ln01, one.u)),
                                                ncol = ln01))
	            if(log) log(res) else res
	        }else{ # explicit
                    alpha <- 1/theta
                    psiI.alpha <- psiI.^alpha
                    sum. <- psiDabsGpoly(psiI.alpha, alpha, d)
                    l.u. <- rowSums(l.u)
                    ll.u. <- rowSums(ll.u)
                    res[n01] <- -psiI.alpha-d*log(psiI.)+ log(sum.) +
                        d*log(theta)+ (theta-1)*ll.u. + l.u.
	            if(log) res else exp(res)
	        }
	    } else res # all NA
	},
        ## parameter interval
        paraInterval = interval("[1,Inf)"),
        ## parameter subinterval (meant to be for *robust* optimization, root-finding etc.)
        paraSubInterval = num2interval(c(1, 50)), # 50 corresponds to tau = 0.98
        ## nesting constraint
        nestConstr = function(theta0,theta1) {
            copGumbel@paraConstr(theta0) &&
            copGumbel@paraConstr(theta1) && theta1 >= theta0
        },
        ## V0 with density dV0 and V01 with density dV01 corresponding to
        ## LS^{-1}[exp(-V_0psi_0^{-1}(psi_1(t)))]
        V0 = function(n,theta) {
            if(theta == 1) {
                ## Sample from S(1,1,0,1;1)
                ## with Laplace-Stieltjes transform exp(-t)
                rep.int(1., n)
            } else {
                alpha <- 1/theta
                ## Sample from S(alpha,1,(cos(alpha*pi/2))^(1/alpha),0;1)
                ## with Laplace-Stieltjes transform exp(-t^alpha)
                rstable1(n, alpha, beta=1,
                         gamma = cos(alpha*pi/2)^(1/alpha))
                ## Note: calling sequence:
                ##       rstable1 -> rstable1C (in rstable1.R) -> rstable1C (in rstable1.R)
                ##       -> rstable_c (in retstable.c) -> rstable_vec (in retstable.c)
                ##       -> rstable0 (in retstable.c)
            }
        },
        dV0 = function(x,theta,log = FALSE) copGumbel@dV01(x,1,1,theta,log),
        V01 = function(V0,theta0,theta1) {
            alpha <- theta0/theta1
            if(alpha == 1) {
                ## Sample from S(1,1,0,V0;1)
                ## with Laplace-Stieltjes transform exp(-V0*t)
                V0
            } else {
                rstable1(length(V0), alpha, beta=1,
                         gamma = (cos(alpha*pi/2)*V0)^(1/alpha))
                ## Sample from S(alpha,1,(cos(alpha*pi/2)V0)^(1/alpha),0;1)
                ## with Laplace-Stieltjes transform exp(-V0*t^alpha)
            }
        },
        dV01 = function(x,V0,theta0,theta1,log = FALSE){
            stopifnot(length(V0) == 1 || length(x) == length(V0))
            alpha <- theta0/theta1
            beta <- 1
            gamma <- (cos(pi/2*alpha)*V0)^(1/alpha)
            delta <- V0*(alpha == 1)
            mapply(dstable,x,alpha = alpha,beta = beta,gamma = gamma,
                   delta = delta, pm = 1, log = log)
        },
        ## Kendall's tau
        tau = function(theta) { (theta-1)/theta },
        tauInv = function(tau) { 1/(1-tau) },
        ## lower tail dependence coefficient lambda_l
        lambdaL = function(theta) { 0*theta },
        lambdaLInv = function(lambda) {
            if(any(lambda != 0))
                stop("Any parameter for a Gumbel copula gives lambdaL = 0")
            NA * lambda
        },
        ## upper tail dependence coefficient lambda_u
        lambdaU = function(theta) { 2 - 2^(1/theta) },
        lambdaUInv = function(lambda) { 1/log2(2-lambda) }
        ) # {copGumbel}

### ==== Joe, see Nelsen (2007) p. 116, # 6 ====================================

##' Joe object
copJoe <-
    new("acopula", name = "Joe",
        ## generator
        psi = function(t,theta) {
            1 - (-expm1(0-t))^(1/theta)
            ## == 1 - (1-exp(-t))^(1/theta)
        },
        psiInv = function(t,theta) { -log1p(-(1-t)^theta) },
        ## absolute value of generator derivatives
        psiDabs = function(t, theta, degree = 1, n.MC, log = FALSE) {
            if(!(missing(n.MC) || is.null(n.MC))){
                psiDabsMC(t,"Joe",theta,degree,n.MC,log)
            }else{
                if(theta == 1) return(if(log) -t else exp(-t)) # independence
                res <- numeric(n <- length(t))
                res[is0 <- t == 0] <- Inf
                res[isInf <- is.infinite(t)] <- if(log)-Inf else 0
                if(any(n0Inf <- !(is0 | isInf))) {
                    alpha <- 1/theta
                    mt <- -t[n0Inf]
                    let <- log( e.t <- -expm1(mt) )
                    ## let = log(e.t) = log(-expm1(-t.)) = log(-(exp(-t.) - 1)) = log(1 - exp(-t.))
		    sum. <- psiDabsJpoly(mt - let, alpha, degree, log=log)
                    res[n0Inf] <-
                        if(log) log(alpha) + alpha*let + sum.
                        else alpha * e.t^alpha *sum.
                }
		res
            }
        },
        ## derivatives of the generator inverse
        psiInvD1abs = function(t, theta, log = FALSE){
            if(log){
                log(theta)+(theta-1)*log1p(-t)-log1p(-(1-t)^theta)
            }else{
                theta/((1-t)^(1-theta)-(1-t))
            }
        },
	## density
	dacopula = function(u, theta, n.MC, log = FALSE,
                            method = c("logJpoly","maxScale","Jpoly")) {
	    if(!is.matrix(u)) u <- rbind(u)
	    if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
            ## f() := NaN outside and on the boundary of the unit hypercube
	    res <- rep.int(NaN, n <- nrow(u))
            n01 <- apply(u,1,function(x) all(0 < x, x < 1)) # indices for which density has to be evaluated
	    if(!any(n01)) return(res)

            ## else - there are u's inside the unit hypercube:
            if(theta == 1){ res[n01] <- if(log) 0 else 1; return(res) } # independence

            ## auxiliary results
            u. <- u[n01,, drop=FALSE]
            u.. <- -(1-u.)^theta # -(1-u)^theta for the relevant u[]
            lpu <- log1p(u..)
            l.x <- rowSums(lpu) # rowSums(log(1-(1-u)^theta)) = log(x)
            l1_u <- log1p(-u.)
            ## main part
            if(!(missing(n.MC) || is.null(n.MC))){ # Monte Carlo
                V <- copJoe@V0(n.MC, theta)
                l <- d*log(theta*V)
                sum. <- (theta-1)*rowSums(l1_u)
                ln01 <- sum(n01)
                one.u <- function(i) exp(l + (V-1)*l.x[i] + sum.[i])
                res[n01] <- colMeans(matrix(unlist(lapply(1:ln01, one.u)),
                                            ncol = ln01))
                if(log) log(res) else res
            }else{
                alpha <- 1/theta
                l.m.x <- log1p(-apply(1 + u.., 1, prod)) # log(1-x)
                l.xx <- l.x - l.m.x # log(x) - log(1-x)
                switch(match.arg(method),
                       "maxScale" = # explicit - scale via factoring out max
                   {
                       ## (1) ingredients
                       l.Stir <- log(Stirling2.all(d))
                       l.gamma <- lgamma((1:d)-alpha)
                       ## (2) carefully compute the sum
                       ## (2.1) compute the values of f
                       ##' f(k), k in 1:d; vector of length n { = length(l.xx) }:
                       f <- function(k) l.Stir[k]+l.gamma[k]-l.gamma[1]+(k-1)*l.xx
                       f.vals <- cbind(0, do.call(cbind, lapply(2:d,f)))
                       ## (2.2) for each i in 1:n, find the k_i^star that maximizes f.vals[i,]
                       k.star <- apply(f.vals, 1, which.max) # k.star[i] = index of maximum of f.vals[i,]
                       f.max <- unlist(lapply(1:n, function(i) f.vals[i,k.star[i]]))
                       ## (2.3) compute the sum
                       f.vals.wo.max.vec <- unlist(lapply(1:n, function(i) f.vals[i,-k.star[i]]))
                       f.vals.wo.max <- matrix(f.vals.wo.max.vec, nrow = n, byrow = TRUE) # f.vals without maxima
                       sum. <- rowSums(exp(f.vals.wo.max))
                       ## (3) compute the density
                       res[n01] <- (d-1)*log(theta) + (theta-1)*rowSums(l1_u) -
                           (1-alpha)*l.m.x + f.max + log1p(sum.)
                   },
                       "logJpoly" =
                   {
                       lsum <- psiDabsJpoly(l.xx, alpha, d, log = TRUE)
                       res[n01] <- (d-1)*log(theta) + alpha*l.m.x + lsum +
                           rowSums((theta-1)* l1_u - lpu)
                   },
                       "Jpoly" =
                   {
                       ## old approach: ok, but numerical problems for large d (~ 100):
                       sum. <- psiDabsJpoly(l.xx, alpha, d)
                       res[n01] <- (d-1)*log(theta) + alpha*l.m.x + log(sum.) +
                           rowSums((theta-1)* l1_u - lpu)
                   }
                       )
                if(log) res else exp(res)
            }
	},
        ## parameter interval
        paraInterval = interval("[1,Inf)"),
        ## parameter subinterval (meant to be for *robust* optimization, root-finding etc.)
        paraSubInterval = num2interval(c(1, 98)), # 98 corresponds to tau = 0.98
        ## nesting constraint
        nestConstr = function(theta0,theta1) {
            copJoe@paraConstr(theta0) &&
            copJoe@paraConstr(theta1) && theta1 >= theta0
        },
        ## V0 with density dV0 and V01 with density dV01 corresponding to
        ## LS^{-1}[exp(-V_0psi_0^{-1}(psi_1(t)))]
        V0 = function(n,theta) rSibuya(n, 1/theta),
        dV0 = function(x,theta,log = FALSE){
            if(log) lchoose(1/theta,x) else abs(choose(1/theta,x))
        },
        V01 = function(V0, theta0, theta1, approx = 10000) {
            ## approx is the largest number of summands before asymptotics is used
            alpha <- theta0/theta1
            rF01Joe(V0, alpha, approx)
        },
        dV01 = function(x,V0,theta0,theta1,log = FALSE){
            l.x <- length(x)
            l.V0 <- length(V0)
            stopifnot(l.V0 == 1 || l.x == l.V0, all(x >= V0))
            alpha <- theta0/theta1
            res <- numeric(l.x)
            if(l.V0 == 1) V0 <- rep(V0,l.x)
            V0.one <- V0 == 1
            res[V0.one] <- copJoe@dV0(x[V0.one],1/alpha,log) # V0 = 1 (numerically more stable)
            if(any(!V0.one)){
                ## FIXME: general case; numerically critical, see, e.g., dV01(1000,500,3,5) < 0
                one.d.args <- function(x.,V0.){
                    j <- 1:V0. # indices of the summands
                    signs <- (-1)^(j+x.)
                    ## determine the signs of choose(j*alpha,x.) for each component of j
                    to.subtract <- 0:(x.-1)
                    signs.choose <- unlist(lapply(j,function(l){
                        prod(sign(l*alpha-to.subtract))}
                                                  ))
                    signs <- signs*signs.choose
                    binom.coeffs <- exp(lchoose(V0.,j)+lchoose(j*alpha,x.))
                    sum(signs*binom.coeffs)
                }
                sum. <- mapply(one.d.args,x[!V0.one],V0[!V0.one])
                res[!V0.one] <- if(log) log(sum.) else sum.
            }
            res
        },
        ## Kendall's tau
        ## noTerms: even for theta==0, the approximation error is < 10^(-5)
        tau = function(theta, noTerms=446) {
            k <- noTerms:1
            sapply(theta,
                   function(th) {
                       tk2 <- th*k + 2
                       1 - 4*sum(1/(k*tk2*(tk2 - th)))
                       ## ==... (1/(k*(th*k+2)*(th*(k-1)+2)))
                   })
        },
        tauInv = function(tau, tol = .Machine$double.eps^0.25, ...) {
            sapply(tau,function(tau) {
                r <- safeUroot(function(th) copJoe@tau(th) - tau,
                               interval = as.numeric(copJoe@paraSubInterval),
                               Sig = +1, tol=tol, check.conv=TRUE, ...)
                r$root
            })
        },
        ## lower tail dependence coefficient lambda_l
        lambdaL = function(theta) { 0*theta },
        lambdaLInv = function(lambda) {
            if(any(lambda != 0))
                stop("Any parameter for a Joe copula gives lambdaL = 0")
            NA * lambda
        },
        ## upper tail dependence coefficient lambda_u
        lambdaU = function(theta) { 2-2^(1/theta) },
        lambdaUInv = function(lambda) { log(2)/log(2-lambda) }
        ) ## {copJoe}

### ==== Generalized inverse Gaussian family ===================================

##' GIG object
                                        # copGIG <-
                                        #     new("acopula", name = "GIG",
                                        #         ## generator
                                        #         psi = function(t,theta) { # theta[1] (nu) in IR, theta[2] (theta) in (0,infty); note that the case theta[1] != 0 and theta[2] == 0 is excluded since this reduces to psiClayton(t,1/theta[1])
                                        #             (1+t)^(-theta[1]/2)*besselK(theta[2]*sqrt(1+t),theta[1])/besselK(theta[2],theta[1]) # theta[1] = nu, theta[2] = theta
                                        #         },
                                        #         psiInv = function(t,theta,...) { ## FIXME: "validTheta" has to be changed (I guess) to allow for more parameters
                                        #             if(is.null(t)) return(as.numeric(NULL)) # FIXME: only a hack to check if validity proceeds
                                        #             r <- safeUroot(function(t) copGIG@psi(t,theta) - t,
                                        #                            interval = c(0,t^(theta[1]/2)-1), Sig = -1, check.conv=TRUE,...)
                                        #             r$root
                                        # 	},
                                        #         ## absolute value of generator derivatives
                                        #         psiDabs = function(t, theta, degree = 1, n.MC, log = FALSE){
                                        #             if(!(missing(n.MC) || is.null(n.MC))){
                                        #                 copGIG@psiDabsMC(t,"GIG",theta,degree,n.MC,log)
                                        #             }else{
                                        # 		res <- numeric(n <- length(t))
                                        #                 res[is0 <- t == 0] <- Inf
                                        #                 res[isInf <- is.infinite(t)] <- 0
                                        #                 n0Inf <- (1:n)[!(is0 | isInf)]
                                        #                 if(length(n0Inf) > 0){
                                        #                     res[n0Inf] <- d*log(theta/2)-theta[2]*(sqrt(1+t)-1)+
                                        #                         log(besselK(theta[2]*sqrt(1+t),theta[1],expon.scaled=TRUE))-
                                        #                             log(besselK(theta[2],theta[1],expon.scaled=TRUE))-((theta[1]+d)/2)*log1p(t)
                                        #                 }
                                        # 		if(log) res else exp(res)
                                        #             }
                                        #         },
                                        #         ## derivatives of the generator inverse
                                        #         psiInvD1abs = function(t, theta, log = FALSE){
                                        #             res <- -log(copGIG@psiDabs(copGIG@psiInv(t,theta),theta,log = TRUE))
                                        #             if(log) res else exp(res)
                                        #         },
                                        #         ## parameter interval
                                        #         paraInterval = interval("(-Inf,Inf)"), # FIXME: theta[2] has to be in (0,Inf)
                                        #         ## parameter subinterval (meant to be for *robust* optimization, root-finding etc.)
                                        #         paraSubInterval = num2interval(c(-100,100)), # FIXME what to really put in?
                                        #         ## nesting constraint
                                        #         nestConstr = function(theta0,theta1) {
                                        #             stop("It is currently not known if GIG generators can be nested.")
                                        #         },
                                        #         ## V0 with density dV0 and V01 with density dV01 corresponding to
                                        #         ## LS^{-1}[exp(-V_0psi_0^{-1}(psi_1(t)))]
                                        #         V0 = function(n,theta){
                                        #             dens <- udgig(theta[1],1,theta[2]^2) # three args: lambda, psi, chi for the GIG density as on page 497 in McNeil, Frey, Embrechts (2005)
                                        #             gen <- pinvd.new(dens) # approximation of the quantile function by piecewise polynomials
                                        #             ur(gen,n)/2
                                        #         },
                                        #         dV0 = function(x,theta,log = FALSE){
                                        #             dens <- udgig(theta[1],1,theta[2]^2)
                                        #             if(log) log(ud(dens,x)) else ud(dens,x)
                                        #         },
                                        #         V01 = function(V0,theta0,theta1) {
                                        #             stop("It is currently not known if GIG generators can be nested.")
                                        #         },
                                        #         dV01 = function(x,V0,theta0,theta1,log = FALSE){
                                        #             stop("It is currently not known if GIG generators can be nested.")
                                        #         },
                                        #         ## Kendall's tau
                                        #         tau = function(theta,...) {
                                        # 	    integrand <- function(t,theta){
                                        # 		t*besselK(theta[2]*sqrt(1+t),theta[1])^2/(1+t)^(theta[1]+1)
                                        #             }
                                        #             one.th <- function(theta,...){
                                        #                 1-(theta[2]/besselK(theta[2],theta[1]))^2*
                                        #                     integrate(integrand,lower=0,upper=Inf,...)
                                        #             }
                                        #             unlist(lapply(theta,one.th,...))
                                        #         },
                                        #         tauInv = function(tau, tol = .Machine$double.eps^0.25, ...) {
                                        #                                         # sapply(tau,function(tau) {
                                        #                                         #                 r <- safeUroot(function(th) copGIG@tau(th) - tau,
                                        #                                         #                                interval = as.numeric(copGIG@paraSubInterval),
                                        #                                         #                                Sig = +1, tol=tol, check.conv=TRUE, ...)
                                        #                                         #                 r$root
                                        #                                         #             })
                                        #             stop("tauInv() is currently not implemented for copGIG") # FIXME: which parameter to solve for?
                                        #         },
                                        #         ## lower tail dependence coefficient lambda_l
                                        #         lambdaL = function(theta) { 0*theta },
                                        #         lambdaLInv = function(lambda) {
                                        #             if(any(lambda != 0))
                                        #                 stop("Any parameter for a GIG copula gives lambdaL = 0")
                                        #             NA * lambda
                                        #         },
                                        #         ## upper tail dependence coefficient lambda_u
                                        #         lambdaU = function(theta) { 0*theta },
                                        #         lambdaUInv = function(lambda) {
                                        #             if(any(lambda != 0))
                                        #                 stop("Any parameter for a GIG copula gives lambdaU = 0")
                                        #             NA * lambda
                                        #         })

### ==== naming stuff ==========================================================

cNms <- c("copAMH", "copClayton", "copFrank", "copGumbel", "copJoe") #, "copGIG")
## == dput(ls("package:nacopula",pat="^cop"))
nmsC <- unlist(lapply(cNms, function(.)get(.)@name))
sNms <- abbreviate(nmsC, 1)
## keep these {hidden, for now}:
c_shortNames <- structure(sNms, names = nmsC)
c_longNames  <- structure(nmsC, names = sNms)
c_objNames   <- structure(cNms, names = nmsC)
rm(cNms, nmsC, sNms)

