## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
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


### Goodness-of-fit testing for nested Archimedean copulas

### transformations to univariate quantities ###################################

##' Kendall distribution function
##'
##' @title Kendall distribution function
##' @param u evaluation point(s)
##' @param cop acopula with specified parameter
##' @param d dimension
##' @param n.MC if > 0 a Monte Carlo approach is applied with sample size equal
##'	   to n.MC; otherwise the exact formula is used
##' @param log logical indicating whether the logarithm of the Kendall
##'        distribution function is returned
##' @return Kendall distribution function at u
##' @author Marius Hofert
K <- function(u, cop, d, n.MC=0, log=FALSE)
{
    ## checks
    stopifnot(is(cop, "acopula"))
    th <- cop@theta
    ## limiting cases
    n <- length(u)
    res <- numeric(n)
    res[is0 <- u == 0] <- if(log) -Inf else 0
    res[is1 <- u == 1] <- if(log) 0 else 1
    not01 <- seq_len(n)[!(is0 | is1)]
    uN01 <- u[not01]
    ## computations
    if(n.MC > 0) { # Monte Carlo
	stopifnot(is.finite(n.MC))
        if(length(not01)){
            V <- cop@V0(n.MC, th) # vector of length n.MC
            psiI <- cop@psiInv(uN01, th) # vector of length n
            lr <- unlist(lapply(psiI, function(psInv){
                -log(n.MC) + lsum(as.matrix(ppois(d-1, V*psInv, log.p=TRUE)))
                ## Former code: mean(ppois(d-1, V*psInv))
            }))
            res[not01] <- if(log) lr else exp(lr)
        }
    } else { # direct
	if(length(not01))
	    res[not01] <- if(d == 1) { # d == 1
		if(log) log(uN01) else uN01 # K(u) = u
	    } else if(d == 2) { # d == 2
		r <- uN01 + exp( cop@psiInv(uN01, theta=th, log=TRUE) -
                                 cop@psiInvD1abs(uN01, th, log=TRUE) ) # K(u) = u - psi^{-1}(u) / (psi^{-1})'(u)
                if(log) log(r) else r
	    } else { # d >= 3
		j <- seq_len(d-1)
		lpsiI. <- cop@psiInv(uN01, theta=th, log=TRUE)
		lpsiDabs <- do.call(rbind,
				    lapply(j, function(j.)
					   cop@psiDabs(exp(lpsiI.),
						       theta=th,
						       degree=j.,
						       log=TRUE))) # (d-1) x n matrix with n = length(not01) containing log( (-1)^j * psi^{(j)}(psi^{-1}(u)) ) in the j-th row
		lfac.j <- cumsum(log(j)) ## == lfactorial(j)
                lx <- lpsiDabs + j %*% t(lpsiI.) - lfac.j # (d-1) x n matrix
                lx <- rbind(log(uN01), lx) # d x n matrix containing the logarithms of the summands of K
                ls <- lsum(lx) # log(K(u))
                if(log) ls else pmin(1, exp(ls)) # ensure we are in [0,1] {numerical inaccuracy}
		## Former code:
		## K2 <- function(psInv) {
		##    lpsiDabs <- unlist(lapply(j, cop@psiDabs,
		##			      u=psInv, theta=th, log=TRUE))
		##    sum(exp(lpsiDabs + j*log(psInv) - lfac.j))
		## }
		## pmin(1, uN01 + unlist(lapply(psiI[not01], K2)))
		##
		## NB: AMH, Clayton, Frank are numerically not quite monotone near one;
		## --  this does not change that {but maybe slightly *more* accurate}:
		## psiDabs. <- unlist(lapply(j, cop@psiDabs, u = psInv, theta = th,
		##						 log = FALSE))
		##		       sum(psiDabs.*psInv^j/factorial(j))
	    } # else (d >= 3)
    } # if/else method (MC/direct)
    res
}


##' Test statistics for various goodness-of-fit tests of (supposedly) U[0,1]^d
##' distributed vectors of random variates
##'
##' @title Test statistics for U[0,1]^d
##' @param u matrix of (pseudo-/copula-)observations
##' @param method various test statistics. Available are:
##'        "AnChisq": Anderson-Darling test statistic after map to a chi-square
##'                   distribution using the standard normal quantile function
##'        "AnGamma": Anderson-Darling test statistic after map to an Erlang/Gamma
##'                   distribution using the logarithm
##'        "SnB"  : the test statistic S_n^{(B)} in Genest, Remillard, Beaudoin (2009)
##'        "SnC"  : the test statistic S_n^{(C)} in Genest, Remillard, Beaudoin (2009)
##' @return values of the chosen test statistic
##' @author Marius Hofert and Martin Maechler
gofTstat <- function(u, method = c("AnChisq", "AnGamma", "SnB", "SnC"))
{
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    d <- ncol(u)
    n <- nrow(u)
    method <- match.arg(method)
    switch(method,
	   "AnChisq" = ad.test( pchisq(rowSums(qnorm(u)^2), d) )$statistic,
	   "AnGamma" = ad.test( pgamma(rowSums(-log(u)), shape=d) )$statistic,
	   "SnB" =
       { ## S_n(B)
           lu2 <- log1p(-u^2) # n x d matrix of log(1-u_{ij}^2)
           ## Note (modulo rowSums/colSums):
           ## Idea: sum1 = sum(prod(1-u^2)) = sum(exp(sum(lu2)))
           ## = exp(log( sum(exp(rowSums(lu2))) )) = exp(lsum(rowSums(lu2)))
           slu2 <- rowSums(lu2) # vector of length n
	   sum1 <- exp(lsum(matrix(slu2, ncol=1))) # lsum() needs a matrix; result: 1 value
           ## 2) The notation here is similar to Genest, Remillard,
           ## Beaudoin (2009) but we interchange k and j (since j always runs
           ## in 1:d). That being said...
	   lu <- t(log1p(-u)) # t(n x d matrix of log(1-u_{ij})) --> accessing columns
           ln <- log(n)
           ## Idea:
           ##   1/n sum_i sum_k prod_j (1-max(u_{ij},u_{kj}))
           ## = 1/n sum_i sum_k exp( sum_j log(1-max{.,.}) )
           ## = 1/n sum_i sum_k exp( sum_j log(min{1-u_{ij},1-u_{kj}}) )
           ## = 1/n sum_i sum_k exp( sum_j min{ log(1-u_{ij}), log(1-u_{kj}) })
           ## = 1/n sum_i sum_k exp( sum(pmin{ lu[i,], lu[k,]}) )
           ## = 1/n sum_i exp( log(sum_k exp( sum(pmin{ lu[i,], lu[k,]}) )) )
           ## = 1/n sum_i exp( lsum( sum(pmin{ lu[i,], lu[k,]}) ) )
           ## = sum_i exp(-log(n) + lsum( sum(pmin{ lu[i,], lu[k,]}) ))
           ## = sum_i exp(-log(n) + lsum_{over k in 1:n}( sum(pmin{ lu[i,], lu[k,]}) ))
           ## => for each fixed i, (l)apply lsum()
	   sum2mands <- unlist(lapply(1:n, function(i){
	       lu.i <- lu[,i] ## now i is fixed
	       sum.k <- vapply(1:n, function(k)# sum over k (n-dim. vector)
			       sum(pmin(lu.i, lu[,k])), 0.)
	       ls.i <- lsum(matrix(sum.k, ncol=1)) # lsum( sum(pmin(...)) ) for fixed i; 1 value
	       exp(-ln + ls.i)
	   }))
	   n/3^d - sum1/2^(d-1) + sum(sum2mands)
       },
	   "SnC" =
       { ## S_n(C)
	   Dn <- apply(u, 1, function(u.){ # Dn is a vector of length n
	       ## u. is one row. We want to know the number of rows of u
	       ## that are (all) componentwise <= u.
	       mean(apply(t(u) <= u., 2, all)) # TRUE <=> the whole row in u is <= u.
	   })
           Cperp <- apply(u, 1, prod)
	   sum((Dn-Cperp)^2)
       },
	   stop("unsupported method ", method))
}


### multivariate transformations ###############################################

##' Transforms vectors of random variates following the given (nested) Archimedean
##' copula (with specified parameters) to U[0,1]^d vectors of random variates
##' via Rosenblatt's transformation
##'
##' @title Rosenblatt transformation for a (nested) Archimedean copula
##' @param u matrix of (pseudo-/copula-)observations
##' @param cop an outer_nacopula
##' @param m # order up to which Rosenblatt's transform is computed, i.e.,
##'        C(u_j | u_1,...,u_{j-1}), j=2,..,m
##' @param n.MC parameter n.MC for evaluating the derivatives via Monte Carlo
##' @return matrix of supposedly U[0,1]^d realizations
##' @author Marius Hofert
rtrafo <- function(u, cop, m=d, n.MC=0)
{
    d. <- dim(u)
    n <- d.[1]
    d <- d.[2]
    stopifnot(is(cop, "outer_nacopula"), 2 <= m, m <= d)
    if(length(cop@childCops))
        stop("currently, only Archimedean copulas are provided")
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    stopifnot(0 <= u, u <=1)
    cop <- cop@copula
    th <- cop@theta
    stopifnot(cop@paraConstr(th))
    psiI <- cop@psiInv(u, theta=th)
    psiI. <- t(apply(psiI, 1, cumsum))
    ## compute all conditional probabilities
    if(n.MC==0){
        ## Note: C(u_j | u_1,...,u_{j-1}) = \psi^{(j-1)}(\sum_{k=1}^j \psi^{-1}(u_k)) / \psi^{(j-1)}(\sum_{k=1}^{j-1} \psi^{-1}(u_k))
	C.j <- function(j){ # computes C(u_j | u_1,...,u_{j-1}) with the same idea as for cacopula
	    logD <- cop@psiDabs(as.vector(psiI.[,c(j,j-1)]), theta=th,
                                degree=j-1, n.MC=0, log=TRUE)
            exp(logD[1:n]-logD[(n+1):(2*n)])
        }
    }else{ # n.MC > 0
	## draw random variates
	V <- cop@V0(n.MC, th)
        C.j <- function(j){ # computes C(u_j | u_1,...,u_{j-1}) with the same idea as for cacopula
            ## use same idea as default method of psiDabsMC
            ## only difference: only draw V's once
            arg <- c(psiI.[,j], psiI.[,j-1])
            iInf <- is.infinite(arg)
            logD <- numeric(2*n)
            logD[iInf] <- -Inf
            if(any(!iInf)) logD[!iInf] <- lsum(-V %*% t(arg[!iInf]) +
                                               (j-1) * log(V) - log(n.MC))
            exp(logD[1:n]-logD[(n+1):(2*n)])
        }
    }
    cbind(u[,1],sapply(2:m, C.j))
}

##' Transforms vectors of random variates following the given (nested) Archimedean
##' copula (with specified parameters) to U[0,1]^d (or U[0,1]^(d-1)) vectors of
##' random variates via the transformation of Hering and Hofert (2011)
##'
##' @title Transformation of Hering and Hofert (2011)
##' @param u data matrix in [0,1]^d
##' @param cop an outer_nacopula
##' @param include.K boolean indicating whether the last component, K, is also
##'        used (include.K = TRUE)
##' @param n.MC parameter n.MC for K
##' @param inverse logical indicating whether the inverse of htrafo is computed,
##'        that is, the transformation of Wu, Valdez, Sherris (2006).
##' @return matrix of supposedly U[0,1]^d realizations
##' @author Marius Hofert and Martin Maechler
htrafo <- function(u, cop, include.K=TRUE, n.MC=0, inverse=FALSE)
{
    ## checks
    stopifnot(is(cop, "outer_nacopula"))
    if(length(cop@childCops))
        stop("currently, only Archimedean copulas are provided")
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    stopifnot((d <- ncol(u)) >= 2,
	      0 <= u, u <= 1)
    ## trafos
    th <- cop@copula@theta
    if(inverse){ # "simulation trafo" of Wu, Valdez, Sherris (2006)
        stop("TODO Marius")
    } else { # "goodness-of-fit trafo" of Hofert and Hering (2011)

        lpsiI <- cop@copula@psiInv(u, th, log=TRUE) # matrix log(psi^{-1}(u))
        lcumsum <- matrix(unlist(lapply(1:d, function(j)
                                        lsum(t(lpsiI[,1:j, drop=FALSE])))),
                          ncol=d)
        u. <- matrix(unlist(lapply(1:(d-1), function(k) exp(k*(lcumsum[,k]-
                                                               lcumsum[,k+1])) )),
                     ncol=d-1) # transformed components (uniform under H_0)
        if(include.K) u. <- cbind(u., K(cop@copula@psi(exp(lcumsum[,d]), th),
                                        cop=cop@copula, d=d, n.MC=n.MC))
        u.
    }
}


### Gof wrapper ################################################################

##' Conducts a goodness-of-fit test for the given H0 copula cop based on the
##' (copula) data u
##'
##' @title Goodness-of-fit testing for (nested) Archimedean copulas
##' @param u (copula-)data matrix
##' @param cop outer_nacopula with specified H0 parameters
##' @param n.bootstrap number of bootstrap replications
##' @param estimation.method estimation method, see enacopula
##' @param include.K boolean whether K is included in the transformation
##' @param n.MC parameter n.MC for K
##' @param trafo multivariate goodness-of-fit transformation; available are:
##'        "Hering.Hofert" the transformation of Hering, Hofert (2011)
##'        "Rosenblatt" the transformation of Rosenblatt (1952)
##' @param method test statistic for the test of U[0,1]^d; see gofTstat()
##' @param verbose if TRUE, the progress of the bootstrap is displayed
##' @param ... additional arguments to enacopula
##' @return htest object
##' @author Marius Hofert and Martin Maechler
gnacopula <- function(u, cop, n.bootstrap,
		      estimation.method= eval(formals(enacopula)$method),
		      include.K= TRUE, n.MC= 0,
		      trafo= c("Hering.Hofert", "Rosenblatt"),
		      method= eval(formals(gofTstat)$method),
		      verbose= TRUE, ...)
{
    ## setup
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    stopifnot(0 <= u, u <= 1, is(cop, "outer_nacopula"), (d <- ncol(u)) >= 2,
              max(cop@comp) == d, n.bootstrap >= 0, n.MC >= 0)
    if(length(cop@childCops))
	stop("currently, only Archimedean copulas are provided")
    if(n.bootstrap == 0)
	stop("Choose a reasonable number of bootstrap replications or
apply the transformations yourself,  see ?gnacopula.")
    u.name <- deparse(substitute(u))

    ## additional warnings for now
    estimation.method <- match.arg(estimation.method)
    if(estimation.method != "mle"){
	if(estimation.method == "smle") warning("'estimation.method = \"smle\"' may be time-consuming!") else
	warning("Consistency for the chosen estimation.method is not clear. Additionally, numerical problems might appear.")
    }

    ## build multivariate transformation
    trafo <- match.arg(trafo)
    method <- match.arg(method)
    gtrafomulti <-
        switch(trafo,
               "Hering.Hofert" = {
                   function(u, cop) htrafo(u, cop=cop, include.K=include.K, n.MC=n.MC)
               },
               "Rosenblatt" = {
                   function(u, cop) rtrafo(u, cop=cop, n.MC=n.MC)
               },
               stop("invalid 'trafo' argument"))

    ## build test statistic function and 'meth' string describing the method
    meth <- paste0("Bootstrapped (B =", n.bootstrap,") test of ")
    meth2 <- paste0(method,", est.method = ", estimation.method)
    meth <-
	switch(method,
	       "AnChisq" =, ## A_n
	       "AnGamma" = { ## A_n with gamma distribution
		   paste0(meth, "Anderson and Darling (with trafo = ",
				    trafo, " and method = ", meth2, ")")
	       },
	       "SnB" =,## S_n(B) and S_n(C)
	       "SnC" = {
		   paste0(meth, meth2," (with trafo = ", trafo, ")")
	       },
	       stop("wrong 'method' argument"))

    ## main part --- Bootstrapping ------------------

    ## (1) estimate the parameter by the provided estimation method and
    ##	   define the estimated copula
    theta.hat <- enacopula(u, cop, method=estimation.method, ...)
    cop.hat <- onacopulaL(cop@copula@name, list(theta.hat, 1:d)) # copula with theta.hat

    ## (2) transform the data with the copula with estimated parameter
    u.prime <- gtrafomulti(u, cop=cop.hat) # transformed data in the unit hypercube
    T <- gofTstat(u.prime, method=method) # transformed data

    ## (3) conduct the Anderson-Darling test or compute the test statistic (depends on method)

    ## (4) conduct the parametric bootstrap
    theta.hat. <- numeric(n.bootstrap) # vector of estimators
    T. <- vector("list", n.bootstrap) # vector of test.stat() results
    if(verbose) pb <- txtProgressBar(max=n.bootstrap, style = 3) # setup progress bar
    for(k in 1:n.bootstrap) {

	## (4.1) sample from the copula with estimated parameter and build
	##	     the corresponding pseudo-observations
	u. <- pobs(rnacopula(nrow(u), cop.hat))

	## (4.2) estimate the parameter by the provided method and define
	##	     the estimated copula
	theta.hat.[k] <- enacopula(u., cop, method=estimation.method, ...)
	cop.hat. <- onacopulaL(cop@copula@name, list(theta.hat.[k], 1:d))

	## (4.3) transform the data with the copula with estimated parameter
	u.prime. <- gtrafomulti(u., cop=cop.hat.)

	## (4.4) compute the test statistic
        T.[[k]] <- gofTstat(u.prime., method=method)
        if(verbose) setTxtProgressBar(pb, k) # update progress bar
    }
    if(verbose) close(pb) # close progress bar

    ## (5) build and return results
    structure(class = "htest",
	      list(p.value= mean(unlist(T.) > T),
                   statistic = T, data.name = u.name,
		   method=meth, estimator=theta.hat,
		   bootStats = list(estimator=theta.hat., statistic=T.)))
}
