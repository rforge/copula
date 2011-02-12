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

#### Goodness-of-fit testing for nested Archimedean copulas

##' Transforms supposedly U[0,1]^d distributed vectors of random variates to
##' U[0,1]-distributed variates to check uniformity in a one-dimensional setup
##'
##' @title Transformation to a one-dimensional testing setting
##' @param u matrix of random variates to be transformed
##' @param method either "normal" (map to a chi-square distribution) or "log"
##'        (map to an Erlang distribution) 
##' @return the supposedly U[0,1] distributed variates
##' @author Marius Hofert
g01 <- function(u, method=c("normal", "log")){
    stopifnot(all(0 <= u, u <= 1))
    if(!is.matrix(u)) u <- rbind(u)
    d <- ncol(u)
    method <- match.arg(method)
    u. <- switch(method,
                 "log" = { pgamma(rowSums(-log(u)),shape=d) },
                 "normal" = { pchisq(rowSums(qnorm(u)^2),d) },
                 stop("wrong choice of method"))
    if(any(is.na(u.))) stop("missing values in u. -- cannot use ad.test()")
    ## check U[0,1] of u.
    ad.test(u.)
}

##' Kendall distribution function
##'
##' @title Kendall distribution function
##' @param t evaluation point(s)
##' @param cop acopula with specified parameter
##' @param d dimension
##' @param n.MC if > 0 a Monte Carlo approach is applied with sample size equal 
##'        to n.MC; otherwise the exact formula is used
##' @return Kendall distribution function at t
##' @author Marius Hofert
K <- function(t, cop, d, n.MC=0)
{
    stopifnot(is(cop, "acopula"))
    psiI <- cop@psiInv(t, th <- cop@theta)
    n <- length(t)
    if(n.MC > 0){
	stopifnot(is.numeric(n.MC), is.finite(n.MC))
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
	    } else {
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
		## K2 <- function(psInv){
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
##' @param x data matrix
##' @param cop an outer_nacopula
##' @param do.pseudo boolean indicating whether to compute the pseudo-observations
##' @param include.K boolean indicating whether the last component, K, is also used or not
##' @param n.MC parameter n.MC for K
##' @return matrix of supposedly U[0,1]^d realizations
##' @author Marius Hofert and Martin Maechler
gnacopulatrafo <- function(x, cop, do.pseudo=FALSE, include.K = ncol(x) <= 5,
                           n.MC = if(ncol(x) <= 10) 0 else 10000){
    stopifnot(is(cop, "outer_nacopula"))
    if(length(cop@childCops))
        stop("currently, only Archimedean copulas are provided")
    if(!is.matrix(x)) x <- rbind(x)
    stopifnot((d <- ncol(x)) >= 2) 
    if(do.pseudo) x <- pobs(x)
    stopifnot(all(0 <= x, x <= 1))
    ## trafo
    th <- cop@copula@theta
    psiI <- cop@copula@psiInv(x, th) # matrix psi^{-1}(x)
    cumsum.psiI <- apply(psiI, 1, cumsum) # rowwise cumulative sums; caution: output is transposed
    u. <- matrix(unlist(lapply(1:(d-1), function(k) (cumsum.psiI[k,]/cumsum.psiI[k+1,])^k)), 
                 ncol=d-1) # transformed components (uniformly under H_0)
    if(include.K) u. <- cbind(u., K(cop@copula@psi(cumsum.psiI[d,], th), 
                                    cop=cop@copula, d=d, n.MC=n.MC))
    u.
}

##' Conducts a goodness-of-fit test for the given H0 copula cop based on the
##' data x
##'
##' @title Goodness-of-fit testing for nested Archimedean copulas
##' @param x data matrix
##' @param cop outer_nacopula with specified H0 parameters
##' @param n.bootstrap number of bootstrap replications
##' @param estimation.method estimation method, see enacopula
##' @param include.K boolean whether K is included in the transformation
##' @param n.MC parameter n.MC for K
##' @param method for g01, see there
##' @param do.pseudo boolean indicating whether to compute the pseudo-observations
##'        of the given data
##' @param do.pseudo.sim boolean indicating whether to compute the pseudo-observations
##'        of the simulated data
##' @param ... additional arguments to enacopula
##' @param verbose if TRUE, the progress of the bootstrap is displayed
##' @return p-value
##' @author Marius Hofert and Martin Maechler
##' FIXME: in case a bootstrap is used, it would be good to return a list also giving the estimated parameter
gnacopula <- function(x, cop, n.bootstrap=0, 
                      estimation.method=eval(formals(enacopula)$method),
                      include.K = ncol(x)<=5, n.MC = if(ncol(x) <= 10) 0 else 10000, 
                      method=c("normal", "log"), do.pseudo=TRUE, 
                      do.pseudo.sim=TRUE, verbose=TRUE, ...)
{
    stopifnot(is(cop, "outer_nacopula"))
    if(length(cop@childCops))
        stop("currently, only Archimedean copulas are provided")
    if(!is.matrix(x)) x <- rbind(x)
    stopifnot((d <- ncol(x)) >= 2) 
    if(do.pseudo) x <- pobs(x)
    stopifnot(all(0 <= x, x <= 1))

    ## main part
    if(n.bootstrap > 0){

        n <- nrow(x)
        ## (1) estimate the parameter by the provided method
        theta.hat <- enacopula(x, cop, method=estimation.method, n.MC=n.MC, 
                               do.pseudo=FALSE, ...) # since pseudo-observations were built earlier
        cop.hat <- onacopulaL(cop@copula@name, list(theta.hat, 1:d)) # copula with theta.hat
        ## (2) transform the data with the copula with estimated parameter and compute
        ##     the test result according to the given method
        test <- g01(gnacopulatrafo(x, cop.hat, do.pseudo=FALSE, include.K=include.K,
                                   n.MC=n.MC), method=method)$statistic # compute the value of the test statistic 
        ## (3) conduct the parametric bootstrap
        test.vec <- numeric(n.bootstrap) # vector of test statistics
        for(k in 1:n.bootstrap){
            u <- rnacopula(n, cop.hat) # sample the copula with estimated parameter
            if(do.pseudo.sim) u <- pobs(u) # compute pseudo-observations if necessary (do *not* compute them in the next line)
            theta.hat.k <- enacopula(u, cop, method=estimation.method, n.MC=n.MC, 
                                     do.pseudo=FALSE, ...) # estimate the copula parameter:
            cop.k <- onacopulaL(cop@copula@name, list(theta.hat.k, 1:d)) # copula with theta.hat.k
            test.vec[k] <- g01(gnacopulatrafo(u, cop.k, do.pseudo=FALSE,
                                              include.K=include.K, n.MC=n.MC), 
                               method=method)$statistic # compute the value of the test statistic
            ## progress output
            if(verbose && k%%10 == 0)
                cat(sprintf("bootstrap progress: %4.1f%%\n", k/n.bootstrap*100))
        }
        ## (4) estimate p-value -- FIXME: return an "htest" result
        mean(test.vec > test)

    }else{ # no bootstrap -- but still a test (!)
        g01(gnacopulatrafo(x, cop, do.pseudo=FALSE, include.K=include.K, n.MC=n.MC),
            method=method)$p.value
    }
}
