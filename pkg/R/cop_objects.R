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
        ## parameter interval
        paraInterval = interval("[0,1)"),
	## absolute value of generator derivatives
	psiDabs = function(t, theta, degree=1, n.MC=0, log=FALSE){
	    if(n.MC > 0){
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
	dacopula = function(u, theta, n.MC=0, log=FALSE){
	    if(!is.matrix(u)) u <- rbind(u)
            if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
            ## f() := NaN outside and on the boundary of the unit hypercube
	    res <- rep.int(NaN, n <- nrow(u))
            n01 <- apply(u,1,function(x) all(0 < x, x < 1)) # indices for which density has to be evaluated
            if(!any(n01)) return(res)
            if(theta == 0){ res[n01] <- if(log) 0 else 1; return(res) } # independence
            ## auxiliary results
            u. <- u[n01,, drop=FALSE]
            u.. <- -theta*(1-u.)
            sum. <- rowSums(log(u.))
            sum.. <- rowSums(log1p(u..))
            ## main part
            if(n.MC > 0){ # Monte Carlo
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
        },
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
        tauInv = function(tau, tol=.Machine$double.eps^0.25, ...) {
            if(any(tau > 1/3))
                stop("Impossible for AMH copula to attain a Kendall's tau larger than 1/3")
            sapply(tau,function(tau) {
                r <- safeUroot(function(th) tauAMH(th) - tau,
                               interval = c(0, 1-1e-12),
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
        ## parameter interval
        paraInterval = interval("(0,Inf)"),
        ## absolute value of generator derivatives
        psiDabs = function(t, theta, degree=1, n.MC=0, log=FALSE){
            if(n.MC > 0){
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
	dacopula = function(u, theta, n.MC=0, log=FALSE){
	    if(!is.matrix(u)) u <- rbind(u)
	    if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
            ## f() := NaN outside and on the boundary of the unit hypercube
	    res <- rep.int(NaN, n <- nrow(u))
            n01 <- apply(u,1,function(x) all(0 < x, x < 1)) # indices for which density has to be evaluated
            if(!any(n01)) return(res)
            ## auxiliary results
            u. <- u[n01,, drop=FALSE]
            l.u <- rowSums(-log(u.))
            ## main part
            if(n.MC > 0){ # Monte Carlo
	        l.u.mat <- matrix(rep(l.u, n.MC), nrow=n.MC, byrow=TRUE)
                V <- copClayton@V0(n.MC, theta)
                l <- d*log(theta*V)
                theta. <- 1 + theta
                psiI.sum <- rowSums(copClayton@psiInv(u., theta))
		## stably compute log(colMeans(exp(B))) 
                B <- l + theta.*l.u.mat - outer(V,psiI.sum) # matrix of exponents; dimension n.MC x n ["V x u"]
                max.B <- apply(B, 2, max) 
                res[n01] <- max.B + log(colMeans(exp(B - rep(max.B, each=n.MC))))
            }else{ # explicit
                alpha <- 1/theta
                d.a <- d + alpha
                arg <- rowSums(u.^(-theta)-1)
                res[n01] <- lgamma(d.a)-lgamma(alpha)+ d*log(theta) +
                    (1+theta)*l.u - (d.a)*log1p(arg)
            }
            if(log) res else exp(res)
	},
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
            gamma <- (cos(pi/2*alpha)*V0)^(1/alpha)
            delta <- V0*(alpha == 1)
	    if(FALSE) ## new dstable() is vectorized in (x, gamma, delta) [but not the others]
                dst <- dstable(x, alpha=alpha, beta = 1, gamma=gamma, delta=delta, pm = 1, log=log)
	    else ## old dstable() needs mapply(.)
                dst <- mapply(dstable,x, alpha=alpha, beta = 1, gamma=gamma,
                              delta=delta, pm = 1, log=log)
	    if(log) V0-x + dst else exp(V0-x) * dst
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
        ## parameter interval
        paraInterval = interval("(0,Inf)"),
        ## absolute value of generator derivatives
        psiDabs = function(t, theta, degree=1, n.MC=0, log=FALSE){
            if(n.MC > 0){
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
	dacopula = function(u, theta, n.MC=0, log=FALSE){
	    if(!is.matrix(u)) u <- rbind(u)
	    if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
            ## f() := NaN outside and on the boundary of the unit hypercube
	    res <- rep.int(NaN, n <- nrow(u))
	    n01 <- apply(u,1,function(x) all(0 < x, x < 1)) # indices for which density has to be evaluated
	    if(!any(n01)) return(res)
            ## auxiliary results
            u. <- u[n01,, drop=FALSE]
            u.sum <- rowSums(u.)
            lp <- log1p(-exp(-theta)) # log(p), p = 1-exp(-theta)
            lpu <- log1p(-exp(-theta*u.)) # log(1 - exp(-theta * u))
            lu <- rowSums(lpu)
            ## main part
            if(n.MC > 0){ # Monte Carlo
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
	},
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
            ##        lgammacor gives underflow warnings for rej < 1
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
            res <- numeric(l)
            res[isN <- theta == 0] <- 0 # limiting case
            res[na <- is.na(theta)] <- NA
            res[!(na | isN)] <- 1 + 4*(debye_1(theta[!na]) - 1)/theta[!na]
            res
        },
        tauInv = function(tau, tol = .Machine$double.eps^0.25, ...){
	    res <- tau
	    res[isN <- res == 0] <- 0 # limiting case	    
            res[!isN] <- sapply(res[!isN], function(tau) {
                r <- safeUroot(function(th) copFrank@tau(th) - tau,
                               interval = c(1e-12, 198),
                               Sig = +1, tol=tol,
                               check.conv=TRUE, ...)
                r$root
            })
            res
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
        ## parameter interval
        paraInterval = interval("[1,Inf)"),
        ## absolute value of generator derivatives
        psiDabs = function(t, theta, degree=1, n.MC=0, 
	method=eval(formals(polyG)$method), log = FALSE){
            if(n.MC > 0){
                psiDabsMC(t,"Gumbel",theta,degree,n.MC,log)
            }else{
                if(theta == 1) return(if(log) -t else exp(-t)) # independence
                res <- numeric(n <- length(t))
                res[is0 <- t == 0] <- Inf
                res[isInf <- is.infinite(t)] <- -Inf
                if(any(n0Inf <- !(is0 | isInf))){
	            t. <- t[n0Inf]
                    alpha <- 1/theta
                    lt <- log(t.)
		    res <- -degree*lt -t^alpha + polyG(alpha*lt, alpha=alpha, d=degree, 
                                                       method=method, log=TRUE)
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
	dacopula = function(u, theta, n.MC=0, method=eval(formals(polyG)$method),
        log = FALSE){
	    if(!is.matrix(u)) u <- rbind(u)
	    if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
            ## f() := NaN outside and on the boundary of the unit hypercube
	    res <- rep.int(NaN, n <- nrow(u))
	    n01 <- apply(u,1,function(x) all(0 < x, x < 1)) # indices for which density has to be evaluated
	    if(!any(n01)) return(res)
            if(theta == 1){ res[n01] <- if(log) 0 else 1; return(res) } # independence
            ## auxiliary results
            u. <- u[n01,, drop=FALSE]
            mlu <- -log(u.) # -log(u)
            lmlu <- log(mlu) # log(-log(u))
            ## main part
            if(n.MC > 0){ # Monte Carlo
                psiI. <- rowSums(copGumbel@psiInv(u., theta))
                V <- copGumbel@V0(n.MC, theta)
                l <- d*log(theta*V)
                sum. <- rowSums((theta-1)*lmlu + mlu)
                sum.mat <- matrix(rep(sum., n.MC), nrow=n.MC, byrow=TRUE)
                ## stably compute log(colMeans(exp(B))) 
                B <- l - outer(V, psiI.) + sum.mat # matrix of exponents; dimension n.MC x n ["V x u"]
		max.B <- apply(B, 2, max) 
		res[n01] <- max.B + log(colMeans(exp(B - rep(max.B, each=n.MC))))
                if(log) res else exp(res)
            }else{ # explicit
                alpha <- 1/theta
                ## compute lx = alpha*log(sum(psiInv(u., theta)))
                lx <- alpha*log(rowSums(copGumbel@psiInv(u., theta)))
                ## ==== former version [start] (numerically slightly more stable but slow) ====
		## im <- apply(u., 1, which.max)
		## mat.ind <- cbind(seq_len(n), im) # indices that pick out maxima from u.
                ## mlum <- mlu[mat.ind] # -log(u_max)
                ## mlum.mat <- matrix(rep(mlum, d), ncol = d)
                ## lx <- lmlu[mat.ind] + alpha*log(rowSums((mlu/mlum.mat)^theta)) # alpha*log(sum(psiInv(u, theta)))
                ## ==== former version [end] ====
                ## compute sum
                lsum <- polyG(lx, alpha, d, method=method, log=TRUE)-d*lx/alpha
                ## the rest
                cop.val <- pnacopula(onacopulaL("Gumbel", list(theta, 1:d)), u.)
                res[n01] <- d*log(theta) + rowSums((theta-1)*lmlu + mlu) + lsum
                res[n01] <- if(log) log(cop.val) + res[n01] else cop.val * exp(res[n01])
                res
            }
	},
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
            gamma <- (cos(pi/2*alpha)*V0)^(1/alpha)
            delta <- V0*(alpha == 1)
	    if(FALSE) ## new dstable() is vectorized in (x, gamma, delta) [but not the others]
                dstable(x, alpha=alpha, beta = 1, gamma=gamma, delta=delta, pm = 1, log=log)
	    else ## old dstable() needs mapply(.)
                mapply(dstable,x, alpha=alpha, beta = 1, gamma=gamma, delta=delta, pm = 1, log=log)
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
        ## parameter interval
        paraInterval = interval("[1,Inf)"),
        ## absolute value of generator derivatives
        psiDabs = function(t, theta, degree = 1, n.MC=0, 
	method=eval(formals(polyJ)$method), log = FALSE) {
            if(n.MC > 0){
                psiDabsMC(t, "Joe", theta, degree, n.MC, log)
            }else{
                if(theta == 1) return(if(log) -t else exp(-t)) # independence
                res <- numeric(n <- length(t))
                res[is0 <- t == 0] <- Inf
                res[isInf <- is.infinite(t)] <- -Inf
                if(any(n0Inf <- !(is0 | isInf))) {
	            alpha <- 1/theta
                    mt <- -t[n0Inf] # -t
                    l1mt <- log(-expm1(mt)) # log(1-exp(-t))
		    sum. <- polyJ(mt-l1mt, alpha, degree, method=method, log=log)
                    res[n0Inf] <- -log(theta) + mt - (1-alpha)*l1mt + sum.
                }
		if(log) res else exp(res)
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
	dacopula = function(u, theta, n.MC=0, method=eval(formals(polyJ)$method),
        log = FALSE){
	    if(!is.matrix(u)) u <- rbind(u)
	    if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
	    ## f() := NaN outside and on the boundary of the unit hypercube
	    res <- rep.int(NaN, n <- nrow(u))
	    n01 <- apply(u,1,function(x) all(0 < x, x < 1)) # indices for which density has to be evaluated
	    if(!any(n01)) return(res)
	    if(theta == 1){ res[n01] <- if(log) 0 else 1; return(res) } # independence
	    ## auxiliary results
	    u. <- u[n01,, drop=FALSE]
	    l1_u <- log1p(-u.) # log(1-u)
	    u.. <- (1-u.)^theta # (1-u)^theta
	    lpu <- log1p(-u..) # log(1-(1-u)^theta)
	    lh <- rowSums(lpu) # rowSums(log(1-(1-u)^theta)) = log(x)
	    ## main part
	    if(n.MC > 0){ # Monte Carlo
	        V <- copJoe@V0(n.MC, theta)
	        l <- d*log(theta*V)
	        sum. <- (theta-1)*rowSums(l1_u)
	        sum.mat <- matrix(rep(sum., n.MC), nrow=n.MC, byrow=TRUE)
	        ## stably compute log(colMeans(exp(B))) 
		B <- l + outer(V-1, lh) + sum.mat # matrix of exponents; dimension n.MC x n ["V x u"]
		max.B <- apply(B, 2, max) 
		res[n01] <- max.B + log(colMeans(exp(B - rep(max.B, each=n.MC))))
	    }else{
	        alpha <- 1/theta
	        h <- apply(1 - u.., 1, prod) # h(u) = \prod_{j=1}^d (1-(1-u_j)^\theta)
	        l1_h <- log1p(-h) # log(1-h(u))
	        lh_l1_h <- lh - l1_h # log(h(u)/(1-h(u)))
		res[n01] <- (d-1)*log(theta) + (theta-1)*rowSums(l1_u) - 
                    (1-alpha)*log1p(-h) + polyJ(lh_l1_h, alpha, d, method=method, 
                                                log=TRUE)
            }
            if(log) res else exp(res)
        },
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
                               interval = c(1, 98),
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

