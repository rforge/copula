## Copyright (C) 2010--2011  Marius Hofert and Martin Maechler
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

#### Goodness-of-fit testing for nested Archimedean copulas

##' Transforms supposedly U[0,1]^d distributed vectors of random variates to
##' U[0,1]-distributed variates (for checking uniformity in a one-dimensional
##' setup)
##'
##' @title Transformation to a one-dimensional test setting
##' @param u matrix of random variates to be transformed
##' @param method either "normal" (map to a chi-square distribution) or "log"
##'	   (map to an Erlang distribution)
##' @return the supposedly U[0,1] distributed variates
##' @author Marius Hofert
g01trafo <- function(u, method = c("normal", "log")) {
    stopifnot(0 <= u, u <= 1)
    if(!is.matrix(u)) u <- rbind(u)
    d <- ncol(u)
    method <- match.arg(method)
    switch(method,
	   "log" = pgamma(rowSums(-log(u)), shape=d),
	   "normal" = pchisq(rowSums(qnorm(u)^2), d),
	   stop("unsupported method ", method))
}

##' Anderson-Darling test for U[0,1] of supposedly U[0,1]^d distributed data
##'
##' @title Anderson-Darling test for U[0,1]
##' @param u matrix of random variates
##' @param method either "normal" (map to a chi-square distribution) or "log"
##'	   (map to an Erlang distribution)
##' @return Anderson-Darling test result
##' @author Marius Hofert
g01 <- function(u, method=c("normal", "log")) {
    u. <- g01trafo(u, method=method)
    if(any(is.na(u.))) stop("missing values in u. -- cannot use ad.test()")
    ad.test(u.) # check U[0,1] for u.
}

##' Kendall distribution function
##'
##' @title Kendall distribution function
##' @param t evaluation point(s)
##' @param cop acopula with specified parameter
##' @param d dimension
##' @param n.MC if > 0 a Monte Carlo approach is applied with sample size equal
##'	   to n.MC; otherwise the exact formula is used
##' @return Kendall distribution function at t
##' @author Marius Hofert
K <- function(t, cop, d, n.MC=0)
{
    stopifnot(is(cop, "acopula"), is.numeric(n.MC),
	      is.numeric(d), length(d) == 1, d == round(d), d >= 1)
    psiI <- cop@psiInv(t, th <- cop@theta)
    n <- length(t)
    if(n.MC > 0) {
	stopifnot(is.finite(n.MC))
	V <- cop@V0(n.MC,th)
	unlist(lapply(psiI, function(psInv) mean(ppois(d-1, V* psInv))))
    } else {
	K. <- numeric(n)
	K.[is0 <- t == 0] <- 0
	K.[is1 <- t == 1] <- 1
	if(length(not01 <- seq_len(n)[!(is0 | is1)]))
	    K.[not01] <- if(d == 1) {
		t[not01]
	    } else if(d == 2) {
		t[not01] + psiI[not01] / cop@psiInvD1abs(t[not01],th)
	    } else { ## d >= 3 :
		j <- seq_len(d-1)
		lpsiDabs <- do.call(rbind,
				    lapply(j, function(j.)
					   cop@psiDabs(psiI[not01],
						       theta=th,
						       degree=j.,
						       log=TRUE))) # (d-1,n)-matrix with n = length(not01)
		lfac.j <- cumsum(log(j)) ## == lfactorial(j)
		r <- colSums(exp(lpsiDabs + j %*% t(log(psiI[not01])) - lfac.j))
		## ensure we are in [0,1] {numerical inaccuracy}
		pmin(1, t[not01] + r)
		## Former code:
		## K2 <- function(psInv) {
		##    lpsiDabs <- unlist(lapply(j, cop@psiDabs,
		##			      t=psInv, theta=th, log=TRUE))
		##    sum(exp(lpsiDabs + j*log(psInv) - lfac.j))
		## }
		## pmin(1, t[not01] + unlist(lapply(psiI[not01], K2)))
		##
		## NB: AMH, Clayton, Frank are numerically not quite monotone near one;
		## --  this does not change that {but maybe slightly *more* accurate}:
		## psiDabs. <- unlist(lapply(j, cop@psiDabs, t = psInv, theta = th,
		##						 log = FALSE))
		##		       sum(psiDabs.*psInv^j/factorial(j))
	    }
	K.
    }
}

##' Transforms vectors of random variates following the given nested Archimedean
##' copula (with specified parameters) to U[0,1]^d vectors of random variates
##'
##' @title Transformation of Hering, Hofert (2011)
##' @param u data matrix
##' @param cop an outer_nacopula
##' @param include.K boolean indicating whether the last component, K, is also used or not
##' @param n.MC parameter n.MC for K
##' @return matrix of supposedly U[0,1]^d realizations
##' @author Marius Hofert and Martin Maechler
gnacopulatrafo <- function(u, cop, include.K = TRUE, n.MC = 0)
{
    stopifnot(is(cop, "outer_nacopula"))
    if(length(cop@childCops))
	stop("currently, only Archimedean copulas are provided")
    if(!is.matrix(u)) u <- rbind(u)
    stopifnot((d <- ncol(u)) >= 2,
	      0 <= u, u <= 1)
    ## trafo
    th <- cop@copula@theta
    psiI <- cop@copula@psiInv(u, th) # matrix psi^{-1}(u)
    cumsum.psiI <- apply(psiI, 1, cumsum) # rowwise cumulative sums; caution: output is transposed
    u. <- matrix(unlist(lapply(1:(d-1), function(k) (cumsum.psiI[k,]/cumsum.psiI[k+1,])^k)),
		 ncol=d-1) # transformed components (uniformly under H_0)
    if(include.K) u. <- cbind(u., K(cop@copula@psi(cumsum.psiI[d,], th),
				    cop=cop@copula, d=d, n.MC=n.MC))
    u.
}

##' Conducts a goodness-of-fit test for the given H0 copula cop based on the
##' data u
##'
##' @title Goodness-of-fit testing for nested Archimedean copulas
##' @param u data matrix
##' @param cop outer_nacopula with specified H0 parameters
##' @param n.bootstrap number of bootstrap replications
##' @param estimation.method estimation method, see enacopula
##' @param include.K boolean whether K is included in the transformation
##' @param n.MC parameter n.MC for K
##' @param method for g01, see there
##' @param ... additional arguments to enacopula
##' @param verbose if TRUE, the progress of the bootstrap is displayed
##' @return if n.bootstrap==0, then the result of the ad.test() call is returned,
##'	        otherwise, a list of results is returned
##' @author Marius Hofert and Martin Maechler
gnacopula <- function(u, cop, n.bootstrap=0,
		      estimation.method=eval(formals(enacopula)$method),
		      include.K = TRUE, n.MC = 0,
		      method = c("normal", "log"), verbose=TRUE, ...)
{
    u.name <- deparse(substitute(u))
    stopifnot(is(cop, "outer_nacopula"),
	      (d <- NCOL(u)) >= 2, 0 <= u, u <= 1)
    if(length(cop@childCops))
	stop("currently, only Archimedean copulas are provided")

    ## main part
    if(n.bootstrap == 0) { ## no bootstrap,
	## transform the data to supposedly U[0,1]-distributed variates, and
	## return AD test:
	g01(gnacopulatrafo(u, cop, include.K=include.K, n.MC=n.MC),
	    method = method)
    } else { # bootstrap
	## (1) estimate the parameter by the provided method and define the
	##     estimated copula
	if(!is.matrix(u)) u <- rbind(u)
	theta.hat <- enacopula(u, cop, method=estimation.method, n.MC=n.MC, ...)
	cop.hat <- onacopulaL(cop@copula@name, list(theta.hat, 1:d)) # copula with theta.hat
	## (2) transform the data with the copula with estimated parameter
	trafo <- gnacopulatrafo(u, cop.hat, include.K=include.K, n.MC=n.MC)
					# transformed data in the unit hypercube
	Y <- g01trafo(trafo, method=method) # transformed data in the unit interval
	## (3) conduct the Anderson-Darling test
	if(any(is.na(Y))) stop("gnacopula: cannot use ad.test() due to missing values")
	AD.test <- ad.test(Y)

	## (4) conduct the parametric bootstrap
	theta.hat. <- numeric(n.bootstrap) # vector of estimators
	AD.test. <- vector("list", n.bootstrap) # vector of ad.test() results
	for(k in 1:n.bootstrap) {
	    ## (4.1) sample from the copula with estimated parameter and build
	    ##	     the corresponding pseudo-observations
	    u. <- pobs(rnacopula(nrow(u), cop.hat))
	    ## (4.2) estimate the parameter by the provided method and define
	    ##	     the estimated copula
	    theta.hat.[k] <- enacopula(u., cop, method=estimation.method, n.MC=n.MC, ...)
	    cop.hat. <- onacopulaL(cop@copula@name, list(theta.hat.[k], 1:d))
	    ## (4.3) transform the data with the copula with estimated parameter
	    trafo. <- gnacopulatrafo(u., cop.hat., include.K=include.K, n.MC=n.MC)
	    Y. <- g01trafo(trafo., method=method)
	    ## (4.4) conduct the Anderson-Darling test
	    AD.test.[[k]] <- ad.test(Y.)
	    ## progress output
	    if(verbose && k %% ceiling(n.bootstrap/100) == 0)
		cat(sprintf("bootstrap progress: %3.0f%%\n", k/n.bootstrap*100))
	}
	## (5) return results
	p.value <- mean(unlist(lapply(AD.test., function(x) x$statistic)) >
		       AD.test$statistic)
	structure(class = "htest",
		  list(p.value=p.value, statistic = theta.hat, data.name = u.name,
		       method = "Bootstrapped Anderson-Darling test of gnacopulatrafo()rmed data",
		       AD.test=AD.test,
		       bootStats = list(theta.hat=theta.hat., AD.test=AD.test.)))
    }
}
