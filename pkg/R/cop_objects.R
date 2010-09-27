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
	## generator derivatives
	psiD = function(t,theta,ord = 1){
            if(ord == 1){
                expmt <- exp(-t)
                -(1-theta)*expmt/(1-theta*expmt)^2
            }else{ 
                ## words of warning:
                ## - do not use for computing the density of an AC 
                ##   (numerically not reasonable due to psiInv)                
                ## - do not use for too large orders
                ## - is only an approximation (k <- 1:Inf would be "correct")
                ## note: one may still use Monte Carlo for large orders
                ##       or compute the first couple of derivatives with a CAS
		l.th <- log(theta)
		l1m.th <- log1p(-theta)
		k <- 1:2000
                l.pk. <- l1m.th + (k-1)*l.th + ord*log(k)
		summands <- function(t) sum(exp(l.pk. - k*t)) # for a single t
		(-1)^ord*unlist(lapply(t,summands))
		## for MC:
		## V0. <- copAMH@V0(100000,theta)
		## l.V0 <- ord*log(V0.)
		## summands <- function(t) mean(exp(l.V0-V0*t))
		## (-1)^ord*unlist(lapply(t,summands))
            }
	},
        ## parameter interval
        paraInterval = interval("[0,1)"),
        ## nesting constraint
        nestConstr = function(theta0,theta1) {
            copAMH@paraConstr(theta0) &&
            copAMH@paraConstr(theta1) && theta1 >= theta0
        },
        ## V0 and V01
        V0 = function(n,theta) rgeom(n, 1-theta) + 1,
        V01 = function(V0,theta0,theta1) {
            rnbinom(length(V0),V0,(1-theta1)/(1-theta0))+V0
        },
	## conditional distribution function C(v|u) of v given u
	cCdf = function(v,u,theta) {
            (1-theta*(1-v))/(v)*((1-theta)/(1-theta*(2-(u+v)+u*v)+
                                            theta^2*(1-u)*(1-v)))^2
	},
	## Kendall distribution function
	K = function(t,theta,d){
            k <- 2000
            pk <- dgeom(k,1-theta)
            one.t <- function(t) sum(pk*ppois(d-1,k*copAMH@psiInv(t,theta)))
            unlist(lapply(t,one.t))
            ## for MC:
            ## V <- copAMH@V0(100000,theta)
            ## K.fun <- function(t) mean(ppois(d-1,V*copAMH@psiInv(t,theta)))
            ## unlist(lapply(t,K.fun))
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
			       interval = c(1e-12, 1-1e-12), Sig = +1,
                               tol = tol, check.conv=TRUE, ...)
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

### ==== Clayton, see Nelsen (2007) p. 116, #1 (slightly simpler form here) ====

copClayton <-
    new("acopula", name = "Clayton",
        ## generator
        psi = function(t,theta) { (1+t)^(-1/theta) },
        psiInv = function(t,theta) { t^(-theta) - 1 },
	## generator derivatives
	psiD = function(t,theta,ord = 1){
            alpha <- 1/theta
            (-1)^ord*prod(alpha+(0:(ord-1)))*(1+t)^(-(alpha+ord))
	},
        ## parameter interval
        paraInterval = interval("(0,Inf)"),
        ## nesting constraint
        nestConstr = function(theta0,theta1) {
            copClayton@paraConstr(theta0) &&
            copClayton@paraConstr(theta1) && theta1 >= theta0
        },
        ## V0 and V01
        V0 = function(n,theta) { rgamma(n, shape = 1/theta) },
        V01 = function(V0,theta0,theta1) { retstable(alpha=theta0/theta1, V0) },
	## conditional distribution function C(v|u) of v given u
	cCdf = function(v,u,theta) {
            exp(-(1+1/theta)*(log(u^(-theta)+v^(-theta)-1)+theta*log(u)))
	},
	## Kendall distribution function
	K = function(t,theta,d){
            if(d==1){
                t
            }else if(d==2){
                t*(1+(1-t^theta)*1/theta)
            }else{ # d >= 3
                k <- 1:(d-1)
                s.fun <- function(j) log(j+1/theta-1)-log(j)
                s.vec <- numeric(d-1)
                s.vec[1] <- -log(theta)
                for(j in 2:(d-1)) s.vec[j] <- s.vec[j-1] + s.fun(j)
                one.t <- function(t) sum(exp(k*log(1-t^theta)+s.vec[k]))
                t(1+unlist(lapply(t,one.t)))	
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
        ## generator derivatives
        psiD = function(t,theta,ord = 1){
            if(ord == 1){
                expmt <- exp(-t)
                p <- 1-exp(-theta)
                (-1/theta)*p*expmt/(1-p*expmt)
            }else{
                ## words of warning:
                ## - do not use for computing the density of an AC 
                ##   (numerically not reasonable due to psiInv)
                ## - do not use for too large orders
                ## - is only an approximation (k <- 1:Inf would be "correct")
                ## note: one may still use Monte Carlo for large orders
                ##       or compute the first couple of derivatives with a CAS
                l.p <- log1p(-exp(-theta))
                l.th <- log(theta)
                k <- 1:2000
                l.pk. <- k*l.p-l.th + (ord-1)*log(k)
                summands <- function(t) sum(exp(l.pk. - k*t)) # for a single t
                (-1)^ord*unlist(lapply(t,summands))
                ## for MC:
                ## V0. <- copFrank@V0(100000,theta)
                ## l.V0 <- ord*log(V0.)
                ## summands <- function(t) mean(exp(l.V0-V0*t))
                ## (-1)^ord*unlist(lapply(t,summands))
            }
        },
        ## parameter interval
        paraInterval = interval("(0,Inf)"),
        ## nesting constraint
        nestConstr = function(theta0,theta1) {
            copFrank@paraConstr(theta0) &&
            copFrank@paraConstr(theta1) && theta1 >= theta0
        },
        ## V0 (algorithm of Kemp (1981)) and V01
        V0 = function(n,theta) { rlog(n,1-exp(-theta)) },
        V01 = function(V0,theta0,theta1) {
            ## FIXME: how to approximate when V0 large? not yet solved
            ## theoretically
            sapply(lapply(V0, rFFrank, theta0=theta0,theta1=theta1), sum)
        },
        ## conditional distribution function C(v|u) of v given u
        cCdf = function(v,u,theta) {
            e.v <- 1-exp(-theta*v)
            e.u <- exp(-v*u)
            p <- 1-exp(-theta)
            e.v*e.u/(p-(1-e.u)*e.v)
        },
        ## Kendall distribution function
	K = function(t,theta,d){
            k <- 2000
            p <- 1-exp(-theta)
            pk <- p^k/(k*(-log1p(-p))) 
            one.t <- function(t) sum(pk*ppois(d-1,k*copFrank@psiInv(t,theta)))
            unlist(lapply(t,one.t))
            ## for MC:
            ## V <- copFrank@V0(100000,theta)
            ## K.fun <- function(t) mean(ppois(d-1,V*copFrank@psiInv(t,theta)))
            ## unlist(lapply(t,K.fun))
	},
        ## Kendall's tau; debye_1() is from package 'gsl' :
        tau = function(theta) 1 + 4*(debye_1(theta) - 1)/theta,
        tauInv = function(tau, tol = .Machine$double.eps^0.25, ...) {
            sapply(tau, function(tau) {
                r <- safeUroot(function(th) copFrank@tau(th) - tau,
                               interval = c(0.001,100), Sig = +1, tol=tol,
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
        ## generator derivatives
        psiD = function(t,theta,ord = 1){
            if(ord == 1){
                alpha <- 1/theta
                copGumbel@psi(t,theta)*(-alpha)*t^(alpha-1)
            }else{
                ## Monte Carlo
                V0. <- copGumbel@V0(100000,theta)
                l.V0 <- ord*log(V0.)
                summands <- function(t) mean(exp(l.V0-V0*t))
                (-1)^ord*unlist(lapply(t,summands))
            }
        },
        ## parameter interval
        paraInterval = interval("[1,Inf)"),
        ## nesting constraint
        nestConstr = function(theta0,theta1) {
            copGumbel@paraConstr(theta0) &&
            copGumbel@paraConstr(theta1) && theta1 >= theta0
        },
        ## V0 and V01
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
            }
        },
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
        ## conditional distribution function C(v|u) of v given u
        cCdf = function(v,u,theta) {
            cop. <- onacopulaL("Gumbel",list(theta,1:2))
            (1+(log(u)/log(v))^theta)^(1/theta-1)*pnacopula(cop.,c(u,v))/u
        },
        ## Kendall distribution function
	K = function(t,theta,d){
            V <- copGumbel@V0(100000,theta)
            K.fun <- function(t) mean(ppois(d-1,V*copGumbel@psiInv(t,theta)))
            unlist(lapply(t,K.fun))
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
        ## generator derivatives
        psiD = function(t,theta,ord = 1){
            alpha <- 1/theta
            if(ord == 1){
                expmt <- exp(-t)
                -alpha*(1-expmt)^(alpha-1)*expmt
            }else{
                ## words of warning:
                ## - do not use for computing the density of an AC 
                ##   (numerically not reasonable due to psiInv)
                ## - do not use for too large orders
                ## - is only an approximation (k <- 1:Inf would be "correct")
                ##   in contrast to AMH and Frank, the series converges slowly!
                ## note: one may still use Monte Carlo for large orders
                ##       or compute the first couple of derivatives with a CAS
                k <- 1:10000
                l.pk. <- abs(lchoose(alpha,k)) + ord*log(k)
                summands <- function(t) sum(exp(l.pk. - k*t)) # for a single t
                (-1)^ord*unlist(lapply(t,summands))
                ## for MC:
                ## V0. <- copJoe@V0(100000,theta)
                ## l.V0 <- ord*log(V0.)
                ## summands <- function(t) mean(exp(l.V0-V0*t))
                ## (-1)^ord*unlist(lapply(t,summands))
            }
        },
        ## parameter interval
        paraInterval = interval("[1,Inf)"),
        ## nesting constraint
        nestConstr = function(theta0,theta1) {
            copJoe@paraConstr(theta0) &&
            copJoe@paraConstr(theta1) && theta1 >= theta0
        },
        ## V0 and V01
        V0 = function(n,theta) rFJoe(n, 1/theta),
        V01 = function(V0,theta0,theta1, approx=100000) {
            alpha <- theta0/theta1
            ia <- (V0 > approx)
            ie <- !ia
            V01 <- V0 # numeric(length(V0))
            V01[ia] <- V0[ia]^(1/alpha) *
                rstable1(sum(ia),alpha, beta=1,
                         gamma = cos(alpha*pi/2)^(1/alpha),
                         delta = as.integer(alpha==1))
            V01[ie] <- sapply(lapply(V0[ie], rFJoe, alpha=alpha),sum)
            V01
        },
        ## conditional distribution function C(v|u) of v given u
        cCdf = function(v,u,theta) {
            u. <- (1-u)^theta
            v. <- (1-v)^theta
            (1-v.)*(1+v.^((1-u.)/u.))^(1/theta-1)
        },
	## Kendall distribution function
	K = function(t,theta,d){
            k <- 10000
            alpha <- 1/theta
            pk <- choose(alpha,k)*(-1)^(k-1)
            one.t <- function(t) sum(pk*ppois(d-1,k*copJoe@psiInv(t,theta)))
            unlist(lapply(t,one.t))
            ## for MC:
            ## V <- copJoe@V0(100000,theta)
            ## K.fun <- function(t) mean(ppois(d-1,V*copJoe@psiInv(t,theta)))
            ## unlist(lapply(t,K.fun))
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
                               interval = c(1.001, 100), Sig = +1, tol=tol,
                               check.conv=TRUE, ...)
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

cNms <- c("copAMH", "copClayton", "copFrank", "copGumbel", "copJoe")
## == dput(ls("package:nacopula",pat="^cop"))
nmsC <- unlist(lapply(cNms, function(.)get(.)@name))
sNms <- abbreviate(nmsC, 1)
## keep these {hidden, for now}:
c_shortNames <- structure(sNms, names = nmsC)
c_longNames  <- structure(nmsC, names = sNms)
c_objNames   <- structure(cNms, names = nmsC)
rm(cNms, nmsC, sNms)

