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
	psiDAbs = function(t,theta,degree = 1,MC,N,log = FALSE){ 
            if(degree == 1){ # exact for degree == 1
                expmt <- exp(-t)
		th.expmt <- theta*expmt
		if(log) log(1-theta)-t-2*log1p(-th.expmt) else (1-theta)*expmt/
                    (1-th.expmt)^2 
            }else{
		if(MC){ # approximation via MC
                    V0. <- copAMH@V0(N,theta)
                    l.V0. <- degree*log(V0.)
                    summands <- function(t) mean(exp(-V0.*t + l.V0.))
		}else{ # approximation via trunction of the series
                    k <- 1:N
                    l.pk <- copAMH@dV0(k,theta,log = TRUE)
                    sum. <- l.pk + degree*log(k)
                    summands <- function(t) sum(exp(-k*t + sum.)) # for a single t
		}
		res <- unlist(lapply(t,summands))
		if(log) log(res) else res
            }
	},
        ## parameter interval
        paraInterval = interval("[0,1)"),
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
	## conditional distribution function C(v|u) of v given u
	cCdf = function(v,u,theta){
            ifelse(v == 0, Inf,
                   (1-theta*(1-v))/v*((1-theta)/(1-theta*(2-(u+v)+u*v)+theta^2*(1-u)*(1-v)))^2)
	},
	## density of the Archimedean copula
	dAc = function(u,theta,MC,N,log = FALSE){ 
            stopifnot(is.matrix(u), all(0 < u, u <= 1))
            if(theta == 0) res <- rep(0,nrow(u))
            if((d <- ncol(u)) == 2){ # d == 2
                u. <- (1-theta*(1-u[,1]))*(1-theta*(1-u[,2]))
                u.. <- theta*u[,1]*u[,2]
                res <- 3*log((1-theta)/(u.-u..))+log(u.+u..)
            }else{ # d > 2
                n <- nrow(u)
                s.psiInv.u <- rowSums(copAMH@psiInv(u,theta))
                u. <- 1-theta*(1-u)
                u.. <- u*u.
                l.u.. <- rowSums(log(u..))
                l.theta. <- log1p(-theta)
                if(MC){ # approximation via MC
                    V0. <- copAMH@V0(N,theta)
                    l.V0. <- d*log(V0.)
                    one.vector.u <- function(j){ # for the jth row of u
                        mean(exp(l.V0.-V0.*s.psiInv.u[j]))
                    }
                    res <- unlist(lapply(1:n,one.vector.u))+d*l.theta.-l.u..
                }else{ # approximation via truncation of the series
                    k <- 1:N
                    l.theta <- log(theta)
                    sum. <- l.theta + rowSums(log(u/u.))
                    d.l.k <- d*log(k)
                    one.vector.u <- function(j){ # for the jth row of u
                        log(sum(exp(d.l.k+k*sum.[j])))
                    }
                    res <- (d+1)*l.theta-l.theta-l.u..+unlist(lapply(1:n,one.vector.u))
                }
            }
            if(log) res else exp(res)
	},
	## Kendall distribution function
	K = function(t, theta, d, MC, N){
            stopifnot(length(theta) == 1, length(d) == 1)
            if(d==1) { # exact for d == 1
                t
            }else if(d==2){ # exact for d == 2
		ifelse(t == 0, 0,
                       t*(1+(1-theta*(1-t))* copAMH@psiInv(t,theta)/(1-theta)))
            }else{ # d >= 3
		if(MC){ # approximation via MC
                    V <- copAMH@V0(N,theta)
                    K.fun <- function(t) mean(ppois(d-1,V*copAMH@psiInv(t,theta)))
                    unlist(lapply(t,K.fun))
		}else{ # approximation via truncation of the series
                    k <- 1:N
                    l.pk <- copAMH@dV0(k,theta,log = TRUE)
                    one.t <- function(t){
                        if(t == 0){
                            0
                        }else if(t == 1){
                            1
                        }else{
                            sum(exp(l.pk+ppois(d-1,k*copAMH@psiInv(t,theta),
                                               log = TRUE)))
                        }
                    }
                    unlist(lapply(t,one.t))
                }
            }
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
        ## absolute value of generator derivatives
        psiDAbs = function(t,theta,degree = 1,log = FALSE){
            alpha <- 1/theta
	    res <- lgamma(degree+alpha)-(degree+alpha)*log1p(t)-lgamma(alpha)
            if(log) res else exp(res)
        },
        ## parameter interval
        paraInterval = interval("(0,Inf)"),
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
        ## conditional distribution function C(v|u) of v given u
        cCdf = function(v,u,theta){
            one.d.args <- function(v,u){
                if(v == 0) if(u == 0) NaN else 0
                exp(-(1+1/theta)*log((u/v)^theta+1-u^theta))
            }
            mapply(one.d.args,u,v)
        },
        ## density of the Archimedean copula
        dAc = function(u,theta,log = FALSE){
            stopifnot(is.matrix(u), all(0 < u, u <= 1))
            d <- ncol(u)
            alpha <- 1/theta
            d.alpha <- d + alpha
            res <- d*log(theta)+(1-theta)*rowSums(log(u))+lgamma(d.alpha)-lgamma(alpha)-
                (d.alpha)*log1p(-d+rowSums(u^(-theta)))
            if(log) res else exp(res)
        },
        ## Kendall distribution function
        K = function(t,theta,d){
            if(d==1){ # exact for d == 1
                t
            }else if(d==2){ # exact for d == 2
                t*(1+(1-t^theta)*1/theta)
            }else{ # d >= 3
                j <- 1:(d-1)
                s <- log(j) + lbeta(j,1/theta)
                one.t <- function(t) sum(exp(j*log1p(-t^theta)-s))
                t*(1+unlist(lapply(t,one.t)))
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
        psiDAbs = function(t,theta,degree = 1,MC,N,log = FALSE){
            if(degree == 1){ # exact for degree == 1
		e.m.th <- exp(-theta)
                p <- 1-e.m.th
		l.p <- log1p(-e.m.th)
                expmt <- exp(-t)
		if(log) -log(theta)+l.p-t-log1p(-p*expmt) else 1/theta*p*expmt/
                    (1-p*expmt)
            }else{
                if(MC){ # approximation via MC
                    V0. <- copFrank@V0(N,theta)
                    l.V0. <- degree*log(V0.)
                    summands <- function(t) mean(exp(-V0.*t + l.V0.))
                }else{ # approximation via trunction of the series
                    k <- 1:N
                    l.pk <- copFrank@dV0(k,theta,log = TRUE)
                    sum. <- l.pk + degree*log(k)
                    summands <- function(t) sum(exp(-k*t + sum.)) # for a single t
                }
                res <- unlist(lapply(t,summands))
                if(log) log(res) else res
            }
        },
        ## parameter interval
        paraInterval = interval("(0,Inf)"),
        ## nesting constraint
        nestConstr = function(theta0,theta1) {
            copFrank@paraConstr(theta0) &&
            copFrank@paraConstr(theta1) && theta1 >= theta0
        },
        ## V0 (algorithm of Kemp (1981)) with density dV0 and V01 with density  
        ## dV01 corresponding to LS^{-1}[exp(-V_0psi_0^{-1}(psi_1(t)))]
        V0 = function(n,theta) { rlog(n,-expm1(-theta)) },
        dV0 = function(x,theta,log = FALSE){
	    if(any(x != (x <- floor(x + 0.5)))) warning("x must be integer; is rounded with floor(x+0.5) otherwise")
            ## FIXME: dgeom, e.g., uses R_D_nonint_check() and R_D_forceint(): is that possible here as well?
            if(log){
                x*log1p(-exp(-theta))-log(x*theta)
            }else{
                p <- 1-exp(-theta)
                p^x/(x*theta)
            }
        },
        V01 = function(V0,theta0,theta1) {
            ## FIXME: how to approximate when V0 large? not yet solved
            ## theoretically
            sapply(lapply(V0, rFFrank, theta0=theta0,theta1=theta1), sum)
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
        ## conditional distribution function C(v|u) of v given u
        cCdf = function(v,u,theta){
            one.d.args <- function(v,u){
                e.v <- -expm1(-theta*v)
                e.u <- -expm1(-theta*u)
		denom <- -expm1(-theta)-e.u*e.v
		if(denom == 0) Inf
                e.v*(1-e.u)/denom
            }
            mapply(one.d.args,u,v)
        },
        ## density of the Archimedean copula
        dAc = function(u,theta,MC,N,log = FALSE){
            stopifnot(is.matrix(u), all(0 < u, u <= 1))
            l.theta <- log(theta)
            if((d <- ncol(u)) == 2){ # d == 2
                p. <- -exp(-theta)
                res <- l.theta+log1p(p.)-theta*(u[,1]+u[,2])-
                    2*log1p(p.-(1-exp(-theta*u[,1]))*(1-exp(-theta*u[,2])))
            }else{ # d > 2
                n <- nrow(u)
                sum.log1p.exp <- rowSums(log1p(-exp(-theta*u)))
                sum.u.theta <- theta*rowSums(u)
                if(MC){ # approximation via MC
                    V0. <- copFrank@V0(N,theta)
                    l.V0. <- d*log(V0.)
                    s.psiInv.u <- rowSums(copFrank@psiInv(u,theta))
                    one.vector.u <- function(j){ # for the jth row of u
                        mean(exp(l.V0.-V0.*s.psiInv.u[j]))
                    }
                    res <- unlist(lapply(1:n,one.vector.u))+d*l.theta-sum.log1p.exp-sum.u.theta
                }else{ # approximation via truncation of the series
                    k <- 1:N
                    d.l.k <- (d-1)*log(k)
                    l.p <- log1p(-exp(-theta))
                    m.l.p <- (1-d)*l.p
                    one.vector.u <- function(j){
                        log(sum(exp(d.l.k+k*(m.l.p+sum.log1p.exp[j]))))
                    }
                    res <- (d-1)*l.theta-sum.u.theta-sum.log1p.exp+unlist(lapply(1:n,one.vector.u))
                }
            }
            if(log) res else exp(res)	
        },
        ## Kendall distribution function
        K = function(t, theta, d, MC, N){
            stopifnot(length(theta) == 1)
            if(d==1){                       # exact for d == 1
                t
            } else {
                ## fast version of
                ##  ifelse(t == 0, 0,
                ##   ifelse(t == 1, 1, { ... }):
                r <- t
                t <- t[n01 <- (t != 0) & (t != 1)]
                r[n01] <-
                    if(d==2) {              # exact for d == 2
                        e.th <- exp(-theta*t)
                        t + copFrank@psiInv(t,theta) * (1-e.th)/(theta*e.th)
                    }
                    else {                  # d >= 3
                        stopifnot(d >= 3)
                        if(MC){             # approximation via MC
                            V <- copFrank@V0(N,theta)
                            K.1.t <- function(t)
                                mean(ppois(d-1, V*copFrank@psiInv(t,theta)))
                        } else {            # approximation via truncation of the series
                            k <- 1:N
                            l.pk <- copFrank@dV0(k,theta,log = TRUE)
                            K.1.t <- function(t)
                                sum(exp(l.pk +
                                        ppois(d-1, k*copFrank@psiInv(t,theta),
                                              log = TRUE)))
                        }
                        unlist(lapply(t, K.1.t))
                    } ## else d >= 3
                r
            } ## else d >= 2
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
        ## absolute value of generator derivatives
        psiDAbs = function(t,theta,degree = 1,MC,N,log = FALSE){
            n <- length(t)
            if(degree == 1){ # exact for degree == 1
                alpha <- 1/theta
                if(log) log(alpha)+(alpha-1)*log(t)-t^alpha else copGumbel@psi(t,theta)*
                    alpha*t^(alpha-1)
            }else{
                res <- rep(0,n) # for those t which are Inf, 0 should be the result
                res[t == 0] <- Inf # for those t which are 0, Inf should be the result
                ind <- (1:n)[t != 0 & t != Inf] # indices of t for which res has to be computed
                if(MC){ # approximation via MC
                    V0. <- copGumbel@V0(N,theta)
                    l.V0. <- degree*log(V0.)
                    summands <- function(t) mean(exp(-V0.*t + l.V0.))
                }else{ # approximation via truncation of the series
                    ## for a single summand:
                    j <- 0:(degree-1)
                    alpha <- 1/theta
                    sign.k <- function(k){ # for a single k
                        s <- sign(j/k-alpha) # degree terms
                        if(any(s == 0)){
                            0
                        }
                        if(sum(s < 0) %% 2 == 0) s <- 1 else s <- -1
                        if(k %% 2) -s else s
                    }
                    ## now look at all summands:
                    k <- 1:N
                    signs <- sign(unlist(lapply(k,sign.k)))
                    ## function for computing sum_{j=0}^{degree-1} log |j/k-alpha| for a single k
                    exp.sum <- function(k) sum(log(abs((0:(degree-1))/k-alpha)))
                    ## compute the summands for a single t
                    factor1 <- alpha*k-degree
                    factor2 <- (degree-1)*log(k)-lfactorial(k-1)+
                        unlist(lapply(k,exp.sum))
                    summands <- function(t) sum(signs*exp(factor1*log(t)+factor2))
                }
                res[ind] <- unlist(lapply(t[ind],summands))
                if(log) log(res) else res
            }
        },
        ## parameter interval
        paraInterval = interval("[1,Inf)"),
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
        ## conditional distribution function C(v|u) of v given u
        cCdf = function(v,u,theta) {
            one.d.args <- function(v,u){
                ## limiting cases
		if(v == 0 && u == 0) NaN
		else if(theta == 1){
		    if(v == 1) 1
		    if(v > 0 && v < 1 && u == 0) v
		}
		else if(v == 1) { if(u == 1) NaN else 0 }
		else if(u == 0 && 0 < v && v < 1) 0
		else {
		    ## main part
		    cop. <- onacopulaL("Gumbel",list(theta,1:2))
		    (1+(log(u)/log(v))^theta)^(1/theta-1)*pnacopula(cop.,c(u,v))/u
		}
            }
            mapply(one.d.args,u,v)
        },
        ## density of the Archimedean copula
        dAc = function(u,theta,MC,N,log = FALSE){
            stopifnot(is.matrix(u), all(0 < u, u < 1))
            if(theta == 1) res <- rep(0,nrow(u))
            u. <- rowSums(copGumbel@psiInv(u,theta))
            alpha <- 1/theta
            if((d <- ncol(u)) == 2){ # d == 2
                res <- log(theta-1+u.^alpha)+(alpha-2)*log(u.)+(theta-1)*
                    log(-log(u[,1])-log(u[,2]))
            }else{ # d > 2
                ## approximation via either truncation of the series or MC
                l.u <- log(u)
                res <- copGumbel@psiDAbs(u.,theta,d,MC=MC,N=N,log=TRUE)+d*
                    log(theta)+(theta-1)*rowSums(log(-l.u))-rowSums(l.u)
            }
            if(log) res else exp(res)
        },
        ## Kendall distribution function
        K = function(t, theta, d, MC, N) {
            stopifnot(length(theta) == 1, length(d) == 1)
            if(d==1){                   # exact for d == 1
                t
            }else if(d==2){             # exact for d == 2
                one.t <- function(t) if(t == 0) 0 else t*(1-log(t)/theta)
                unlist(lapply(t,one.t))
            }else{                      # d >= 3
                if(MC){                 # approximation via MC
                    V <- copGumbel@V0(N,theta)
                    K.fun <- function(t) mean(ppois(d-1,V*copGumbel@psiInv(t,theta)))
                    unlist(lapply(t,K.fun))
                }else{    # approximation via truncation of the series
                    j <- 1:(d-1)
                    l.j <- lfactorial(j)
                    one.t <- function(t){
                        if(t == 0 || t == 1){
                            0
                        }else{
                            pI <- copGumbel@psiInv(t,theta)
                            sum(exp(unlist(lapply(j,copGumbel@psiDAbs,t=pI,
                                                  theta=theta,MC=MC,N=N,log=TRUE))+
                                    j*log(pI)-l.j))
                        }
                    }
                    t + unlist(lapply(t, one.t))
                }
            }
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
        psiDAbs = function(t,theta,degree = 1,MC,N,log = FALSE){
            alpha <- 1/theta
            if(degree == 1){ # exact for degree == 1
                expmt <- exp(-t)
                if(log) log(alpha)+(alpha-1)*log1p(-expmt)-t else alpha*
                    (1-expmt)^(alpha-1)*expmt
            }else{
                if(MC){ # approximation via MC
                    V0. <- copJoe@V0(N,theta)
                    l.V0. <- degree*log(V0.)
                    summands <- function(t) mean(exp(-V0.*t + l.V0.))
                }else{ # approximation via truncation of the series
                    k <- 1:N
                    l.pk <- copJoe@dV0(k,theta,log = TRUE)
                    sum. <- l.pk + degree*log(k)
                    summands <- function(t) sum(exp(-k*t + sum.)) # for a single t
                }
                res <- unlist(lapply(t,summands))
                if(log) log(res) else res
            }
        },
        ## parameter interval
        paraInterval = interval("[1,Inf)"),
        ## nesting constraint
        nestConstr = function(theta0,theta1) {
            copJoe@paraConstr(theta0) &&
            copJoe@paraConstr(theta1) && theta1 >= theta0
        },
        ## V0 with density dV0 and V01 with density dV01 corresponding to 
        ## LS^{-1}[exp(-V_0psi_0^{-1}(psi_1(t)))]
        V0 = function(n,theta) rFJoe(n, 1/theta),
        dV0 = function(x,theta,log = FALSE){
            if(log) lchoose(1/theta,x) else abs(choose(1/theta,x))
        },
        V01 = function(V0,theta0,theta1, approx=10000) {
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
        ## conditional distribution function C(v|u) of v given u
        cCdf = function(v,u,theta){
            one.d.args <- function(u,v){
                u. <- (1-u)^theta
                v. <- (1-v)^theta
                if(theta != 1 && u. == 0){
                    v.1 <- 1-v.
                    if(v. == 1) 0 else v.1
                }
                (1-v.)*(1+v.^((1-u.)/u.))^(1/theta-1)
            }
            mapply(one.d.args,u,v)
        },
        ## density of the Archimedean copula
        dAc = function(u,theta,MC,N,log = FALSE){
            stopifnot(is.matrix(u), all(0 <= u, u < 1))
            alpha <- 1/theta
            if((d <- ncol(u)) == 2){ # d == 2
                u. <- (1-u[,1])^theta
                v. <- (1-u[,2])^theta
                res <- log(theta-(1-u.)*(1-v.))+(1-alpha)*(log(u.)+log(v.))+
                    (alpha-2)*log(u.+v.-u.*v.)
            }else{ # d > 2
                n <- nrow(u)
                d.l.theta <- d*log(theta)
                sum.log.u. <- rowSums(log1p(-u))
                sum.log.u.. <- rowSums(log1p(-(1-u)^theta))
                if(MC){ # approximation via MC
                    V0. <- copJoe@V0(N,theta)
                    l.V0. <- d*log(V0.)
                    s.psiInv.u <- rowSums(copJoe@psiInv(u,theta))
                    one.vector.u <- function(j){ # for the jth row of u
                        mean(exp(l.V0.-V0.*s.psiInv.u[j]))
                    }
                    res <- unlist(lapply(1:n,one.vector.u))+d.l.theta+(theta-1)*
                        sum.log.u.-sum.log.u..
                }else{ # approximation via truncation of the series
                    k <- 1:N
                    l.c <- lchoose(alpha,k)
                    d.l.k <- d*log(k)
                    one.vector.u <- function(j){
                        log(sum(exp(l.c+d.l.k+(k-1)*sum.log.u..[j])))
                    }
                    res <- d.l.theta+(theta-1)*sum.log.u.+unlist(lapply(1:n,one.vector.u))
                }
            }
            if(log) res else exp(res)
        },
        ## Kendall distribution function
        K = function(t,theta,d,MC,N){
            stopifnot(length(theta) == 1, length(d) == 1)
            if(d==1){ # exact for d == 1
                t
            }else if(d==2){ # exact for d == 2
                mt <- 1-t
                one.t <- function(t){
                    if(t == 0){
                        0
                    }else if(t == 1){
                        1
                    }else{
                        mt <- 1-t
                        t+copJoe@psiInv(t,theta)*mt*(mt^(-theta)-1)/theta
                    }
                }
                unlist(lapply(t,one.t))
            }else{ # d >= 3
                if(MC){ # approximation via MC
                    V <- copJoe@V0(N,theta)
                    K.fun <- function(t) mean(ppois(d-1,V*copJoe@psiInv(t,theta)))
                    unlist(lapply(t,K.fun))
                }else{ # approximation via truncation of the series
                    k <- 1:N
                    l.pk <- copJoe@dV0(k,theta,log = TRUE)
                    one.t <- function(t){
                        if(t == 0){
                            0
                        }else if(t == 1){
                            1
                        }else{
                            sum(exp(l.pk+ppois(d-1,k*copJoe@psiInv(t,theta),
                                               log = TRUE)))
                        }
                    }
                    unlist(lapply(t,one.t))
                }
            }
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

