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
##' @param u matrix of random variates to be transformed
##' @param method either "log" (map to an Erlang distribution) or "normal" (map
##' 	   to a chi-square distribution)
##' @return the supposedly U[0,1] distributed variates
##' @author Marius Hofert
g01 <- function(u, method = c("log","normal")){
    stopifnot(all(0 <= u, u <= 1))
    if(is.vector(u)) u <- matrix(u, nrow = 1)
    d <- ncol(u)
    u. <- switch(method,
                 "log" = { pgamma(rowSums(-log(u)),shape=d) },
                 "normal" = { pchisq(rowSums(qnorm(u)^2),d) },
                 stop("wrong choice of method"))
    if(any(is.na(u.))) stop("missing values in u. -- cannot use ad.test()")
    ## check U[0,1] of u.
    ad.test(u.)
}

##' Kendall distribution function
##' @param t evaluation point(s)
##' @param cop acopula with specified parameter
##' @param d dimension
##' @param MC if provided (and not NULL) psiDabs is evaluated via Monte Carlo
##'        with sample size MC
##' @return Kendall distribution function at t
##' @author Marius Hofert
##' FIXME: maybe with outer()
K <- function(t, cop, d, MC)
{
    stopifnot(is(cop, "acopula"))
    psiI <- cop@psiInv(t,th <- cop@theta)
    n <- length(t)
    if(!(missing(MC) || is.null(MC))){
	stopifnot(is.numeric(MC), is.finite(MC), MC > 0)
        V <- cop@V0(MC,th)
        unlist(lapply(psiI, function(psInv) mean(ppois(d-1, V* psInv))))
    } else {
	K. <- numeric(n)
	K.[is0 <- t == 0] <- 0
	K.[is1 <- t == 1] <- 1
	if(length(not01 <- seq_len(n)[!(is0 | is1)]))
	    K.[not01] <- if(d == 1) {
		t[not01]
	    } else if(d == 2) {
		t[not01]+ psiI[not01] / cop@psiInvD1abs(t[not01],th)
	    } else {
		j <- seq_len(d-1)
                lfac.j <- cumsum(log(j)) ## == lfactorial(j)
		K2 <- function(psInv) {
		    lpsiDabs <- unlist(lapply(j, cop@psiDabs,
					      t = psInv, theta = th, log = TRUE))
		    sum(exp(lpsiDabs + j*log(psInv) - lfac.j))
		    ## NB: AMH, Clayton, Frank are numerically not quite monotone near one;
                    ## --  this does not change that {but maybe slightly *more* accurate}:
		    ## psiDabs. <- unlist(lapply(j, cop@psiDabs, t = psInv, theta = th,
		    ##						 log = FALSE))
		    ##		       sum(psiDabs.*psInv^j/factorial(j))
		}
		## ensure we are in [0,1] {numeric inaccuracy}
		pmin(1, t[not01] + unlist(lapply(psiI[not01], K2)))
	    }
	K.
    }
}

##' Transforms vectors of random variates following the given nested Archimedean
##' copula (with specified parameters) to U[0,1]^d vectors of random variates
##' @param x data matrix
##' @param cop an ("outer_") nacopula
##' @param do.pseudo boolean indicating whether to compute the pseudo-observations
##' @param MC if provided K is evaluated via Monte Carlo with sample size MC
##' @return matrix of supposedly U[0,1]^d realizations
##' @author Marius Hofert & Martin Maechler
gnacopulatrafo <- function(x, cop, MC, do.pseudo = FALSE)
{
    stopifnot(is(cop, "outer_nacopula"))
    if(length(cop@childCops))
        stop("currently, only Archimedean copulas are provided")
    if(is.vector(x)) x <- matrix(x, nrow = 1)
    stopifnot((d <- ncol(x)) >= 2)
    u <- if(do.pseudo){
	pobs(x)
    }else{
	stopifnot(all(0 <= x, x <= 1))
	x
    }
    acop <- cop@copula
    th <- acop@theta
    psiI <- acop@psiInv(u, th)
    cumsum.psiI <- t(apply(psiI,1,cumsum))
    u. <- matrix(, nrow = nrow(u), ncol = d) # for the transformed components (U[0,1]^d under H0)
    for(j in seq_len(d-1)) u.[,j] <- (cumsum.psiI[,j]/cumsum.psiI[,j+1])^j
    u.[,d] <- K(acop@psi(cumsum.psiI[,d], th), acop, d = d, MC = MC)
    u.
}

##' Conducts a goodness-of-fit test for the given H0 copula cop based on the
##' data x
##' @param x data matrix
##' @param cop nacopula with specified H0 parameters
##' @param do.pseudo boolean indicating whether to compute the pseudo-observations
##' @param method for Anderson-Darling, see g01
##' @param MC if provided K is evaluated via Monte Carlo with sample size MC
##' @param bootstrap whether or not a bootstrap is applied
##' @param B number of bootstrap replications
##' @param estimation.method estimation method, see enacopula
##' @param verbose if TRUE, the progress of the bootstrap is displayed
##' @return a test (bootstrap or ad.test) test result
##' @author Marius Hofert & Martin Maechler
gnacopula <- function(x, cop, bootstrap = TRUE, B = 1000, method = c("log","normal"),
                      estimation.method = c("mle.tau.mean",
                      "mle.theta.mean","mle.diag","smle","tau.tau.mean",
                      "tau.theta.mean","dmle","beta"), MC, do.pseudo = TRUE,
                      verbose = TRUE)
{
    stopifnot(is(cop, "outer_nacopula"))
    if(is.vector(x)) x <- matrix(x, nrow = 1)
    stopifnot((d <- ncol(x)) >= 2)
    u <- if(do.pseudo){
	pobs(x)
    }else{
	stopifnot(all(0 <= x, x <= 1))
	x
    }

    if(bootstrap){

        ## bootstrap setup
        n <- nrow(u)
        ## although not recommended, bootstrap can be used with do.pseudo == FALSE
        ## in which case the data should come from the copula directly (so at least
        ## it should be between 0 and 1)
        stopifnot(all(0 <= u, u <= 1))
        ## estimate the parameter by the provided method
        theta.hat <- enacopula(u,cop,estimation.method,MC,do.pseudo = do.pseudo)
        ## define the copula with theta.hat
        cop.hat <- onacopulaL(cop@family,list(theta.hat,1:d))
        ## transform the data with the copula with estimated parameter and compute
        ## the AD test result
        AD <- g01(gnacopulatrafo(u, cop.hat, MC, do.pseudo = do.pseudo),
                  method = method)$p.value

        ## conduct the parametric bootstrap
        AD.vec <- numeric(B) # vector of AD test statistics
        for(k in 1:B) {
            ## do work
            u. <- rnacopula(n,cop.hat) # sample the copula with estimated parameter
            if(do.pseudo) u. <- pobs(u.) # compute pseudo-observations if necessary
            ## (and do *not* compute them again below!)

            ## Estimate the copula parameter:
            theta.B <- enacopula(u., cop, estimation.method, MC, do.pseudo=FALSE)
            ## copula with estimated parameter:
            cop.B <- onacopulaL(cop@family,list(theta.B, 1:d))
            AD.vec[k] <- g01(gnacopulatrafo(u., cop.B, MC = MC, do.pseudo=FALSE),
                             method = method)$p.value
            ## progress output
            if(verbose && k%%10 == 0)
                cat(sprintf("bootstrap progress: %4.1f%%\n", k/B * 100))
        }
        ## estimate p-value --- FIXME: return a "htest" result
        mean(AD.vec > AD)

    } else { # no bootstrap --- but still a test (!)
        g01(gnacopulatrafo(u, cop, MC = MC, do.pseudo = FALSE),
            method = method)$p.value
    }
}

