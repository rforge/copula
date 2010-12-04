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
	psiDabs = function(t, theta, degree = 1, MC, log = FALSE){
	    if(!(missing(MC) || is.null(MC))){
		psiDabsMC(t,"AMH",theta,degree,MC,log)
            }else{
		if(theta == 0) if(log) return(-t) else return(exp(-t)) # special case
		## Note: psiDabs(0, ...) is correct
		arg <- theta*exp(-t)
                if(log){
                    log1p(-theta)-log(theta)+log(unlist(lapply(arg,polylog,s = -degree,
                                                               method = "neg")))
                }else{
                    unlist(lapply(arg,polylog,s = -degree,method = "neg"))*
                        (1-theta)/theta
                }
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
        psiDabs = function(t, theta, degree = 1, MC, log = FALSE){
            if(!(missing(MC) || is.null(MC))){
                psiDabsMC(t,"Clayton",theta,degree,MC,log)
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
        psiDabs = function(t, theta, degree = 1, MC, log = FALSE){
            if(!(missing(MC) || is.null(MC))){
                psiDabsMC(t,"Frank",theta,degree,MC,log)
            }else{
                ## Note: psiDabs(0, ...) is correct
                p <- -expm1(-theta)
                arg <- p*exp(-t)
                if(log){
                    -log(theta)+log(unlist(lapply(arg,polylog,s = -(degree-1),
                                                  method = "neg")))
                }else{
                    unlist(lapply(arg,polylog,s = -(degree-1),method = "neg"))/theta
                }
            }
        },
	## derivatives of the generator inverse
	psiInvD1abs = function(t, theta, log = FALSE){
            if(log) log(theta)-log(expm1(theta*t)) else theta/expm1(theta*t)
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
        psiDabs = function(t, theta, degree = 1, MC, log = FALSE){
            if(!(missing(MC) || is.null(MC))){
                psiDabsMC(t,"Gumbel",theta,degree,MC,log)
            }else{
                if(theta == 1) return(if(log) -t else exp(-t)) # special case
                res <- numeric(n <- length(t))
	        res[is0 <- t == 0] <- Inf
	        res[isInf <- is.infinite(t)] <- -Inf
                n0Inf <- (1:n)[!(is0 | isInf)]
		if(length(n0Inf) > 0){ # compute values for indices n0Inf
	            ## compute factors (inner sum)
	            s <- Stirling1.all(degree) # s(d,1), ..., s(d,d)
                    k <- 1:degree
		    S <- lapply(k, Stirling2.all) # S[[m]][n] contains S(m,n), n = 1,...,m
                    a.k <- unlist(lapply(k, function(k.){
                        j <- k.:degree
                        S. <- unlist(lapply(j, function(i) S[[i]][k.])) # extract a column of Stirling2 numbers
                        sum(theta^(-j)*s[j]*S.)
                    }))
                    ## compute outer sum
                    alpha <- 1/theta
                    one.k <- function(k.) (-1)^(degree-k.)*t[n0Inf]^(k.*alpha)*a.k[k.]
		    sum. <- rowSums(as.matrix(sapply(k, one.k))) # FIXME: as.matrix is needed when there is only one value t
		    res[n0Inf] <- -t[n0Inf]^alpha-degree*log(t[n0Inf])+log(sum.)
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
            mapply(dstable,x,alpha = alpha,beta = beta,gamma = gamma,delta = delta,
                   pm = 1,log = log)
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
        psiDabs = function(t, theta, degree = 1, MC, log = FALSE){
            if(!(missing(MC) || is.null(MC))){
                psiDabsMC(t,"Joe",theta,degree,MC,log)
            }else{
		if(theta == 1) return(if(log) -t else exp(-t)) # special case
                res <- numeric(n <- length(t))
	        res[is0 <- t == 0] <- Inf
	        res[isInf <- is.infinite(t)] <- 0
                n0Inf <- (1:n)[!(is0 | isInf)]
		if(length(n0Inf) > 0){ # compute values for indices n0Inf
	            alpha <- 1/theta
                    S <- Stirling2.all(degree)
                    one.k <- function(k) exp(log(S[k])+lgamma(k-alpha)-
                                              lgamma(1-alpha)-k*t[n0Inf]-(k-alpha)*
                                              log1p(-exp(-t[n0Inf])))
                    res[n0Inf] <- alpha*rowSums(as.matrix(sapply(1:degree, one.k))) # FIXME: as for Gumbel, "as.matrix" is needed
                }
	        if(log) log(res) else res
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
                    signs.choose <- unlist(lapply(j,function(l) prod(sign(l*alpha-to.subtract))))
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
                                        #         psiDabs = function(t, theta, degree = 1, MC, log = FALSE){
                                        #             if(!(missing(MC) || is.null(MC))){
                                        #                 copGIG@psiDabsMC(t,"GIG",theta,degree,MC,log)
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

