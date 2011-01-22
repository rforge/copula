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

#### Estimation for nested Archimedean copulas

## ==== initial interval for optimization procedures ===========================

##' Compute an initial interval for optimization/estimation routines (only a
##' heuristic; if this fails, choose your own interval)
##' Note: In contrast to the slots "paraSubInterval", this function also works
##'       for non-robust methods (i.e., methods that break down when the rather
##'	  large intervals given by paraSubInterval are used)
##'
##' @title Compute initial interval for estimation procedures
##' @param u matrix of realizations following a copula
##' @param family Archimedean family
##' @param h for enlarging the tau-interval
##' @return initial interval which can be used for optimization (e.g., for emle)
##' @author Marius Hofert
paraOptInterval <- function(u, family, h=0.15){
    x <- apply(u,1,max)
    theta.hat.G <- log(ncol(u))/(log(length(x))-log(sum(-log(x)))) # direct formula from edmle for Gumbel
    tau.hat.G <- copGumbel@tau(theta.hat.G)
    copFamily <- getAcop(family)
    I <- copFamily@paraSubInterval
    tau.min <- copFamily@tau(I[1]) # smallest admissible lower bound for copFamily
    tau.max <- copFamily@tau(I[2]) # largest admissible upper bound for copFamily
    l <- max(tau.hat.G - h,tau.min) # admissible lower bound for tau
    u <- min(tau.hat.G + h,tau.max) # admissible upper bound for tau
    c(copFamily@tauInv(l),copFamily@tauInv(u))
}

## ==== Blomqvist's beta =======================================================

##' Compute the sample version of Blomqvist's beta,
##' see, e.g., Schmid and Schmidt (2007) "Nonparametric inference on multivariate
##' versions of Blomqvist’s beta and related measures of tail dependence"
##'
##' @title Sample version of Blomqvist's beta
##' @param u matrix of realizations following the copula
##' @param scaling if TRUE then the factors 2^(d-1)/(2^(d-1)-1) and
##'                2^(1-d) in Blomqvist's beta are omitted
##' @return sample version of multivariate Blomqvist beta
##' @author Marius Hofert
beta.hat <- function(u, scaling = FALSE){
    stopifnot(all(0 <= u, u < 1))
    less.u <- u <= 0.5
    prod1 <- apply( less.u, 1, all)
    prod2 <- apply(!less.u, 1, all)
    b <- mean(prod1 + prod2)
    if(scaling) b else {T <- 2^(ncol(u)-1); (T*b - 1)/(T - 1)} 
}

##' Compute the population version of Blomqvist's beta for Archimedean copulas
##'
##' @title Population version of Blomqvist's beta for Archimedean copulas
##' @param cop acopula to be estimated
##' @param theta copula parameter
##' @param d dimension
##' @param scaling if TRUE then the factors 2^(d-1)/(2^(d-1)-1) and
##'                2^(1-d) in Blomqvist's beta are omitted
##' @return population version of multivariate Blomqvist beta
##' @author Marius Hofert & Martin Maechler
beta. <- function(cop, theta, d, scaling = FALSE) {
    stopifnot(is(cop,"acopula"), cop@paraConstr(theta))
    j <- seq_len(d)
    diags <- cop@psi(j*cop@psiInv(0.5, theta), theta) # compute diagonals
    b <- 1 + diags[d] + if(d < 30) sum((-1)^j * choose(d, j) * diags)
    else sum((-1)^j * exp(lchoose(d, j) + log(diags)))
    if(scaling) b else { T <- 2^(d-1); (T*b - 1)/(T - 1)} 
}

##' Method-of-moment-like estimation of nested Archimedean copulas based on a
##' multivariate version of Blomqvist's beta
##'
##' @title Method-of-moment-like parameter estimation of nested Archimedean copulas
##'        based on Blomqvist's beta
##' @param u matrix of realizations following the copula
##' @param cop outer_nacopula to be estimated
##' @param interval bivariate vector denoting the interval where optimization takes
##'        place
##' @param ... additional parameters for safeUroot
##' @return Blomqvist beta estimator; return value of safeUroot (more or less
##'	    equal to the return value of uniroot)
##' @author Marius Hofert
ebeta <- function(u, cop, interval = paraOptInterval(u, cop@copula@name), ...){
    stopifnot(is(cop, "outer_nacopula"))
    if(length(cop@childCops))
        stop("currently, only Archimedean copulas are provided")
    ## Note: We do not need the constants 2^(d-1)/(2^(d-1)-1) and 2^(1-d) here,
    ##	     since we equate the population and sample versions of Blomqvist's
    ##       beta anyway.
    b.hat <- beta.hat(u, scaling = TRUE)
    d <- ncol(u)
    safeUroot(function(theta){beta.(cop@copula, theta, d, scaling=TRUE) - b.hat},
              interval = interval, Sig = +1, check.conv = TRUE, ...)
}

## ==== Kendall's tau ==========================================================

##' Compute pairwise estimators for nested Archimedean copulas based on Kendall's tau 
##'
##' @title Pairwise estimators for nested Archimedean copulas based on Kendall's tau 
##' @param u matrix of realizations following the copula
##' @param cop outer_nacopula to be estimated
##' @param method tau.mean indicates that the average of the sample versions of
##'               Kendall's tau are computed first and then theta is determined;
##'               theta.mean stands for first computing all Kendall's tau
##'               estimators and then returning the mean of these estimators
##' @param ... additional arguments to cor()
##' @return averaged pairwise cor() estimators
##' @author Marius Hofert
etau <- function(u, cop, method = c("tau.mean", "theta.mean"), ...)
{
    stopifnot(is(cop, "outer_nacopula"))
    if(length(cop@childCops))
        stop("currently, only Archimedean copulas are provided")
    method <- match.arg(method)
    tau.hat.mat <- cor(u, method="kendall",...) # matrix of pairwise tau()
    tau.hat <- tau.hat.mat[upper.tri(tau.hat.mat)] # all tau hat's
    switch(method,
           "tau.mean" = {
               cop@copula@tauInv(mean(tau.hat)) # Kendall's tau corresponding to the mean of the tau hat's
           },
           "theta.mean" = {
               mean(cop@copula@tauInv(tau.hat)) # mean of the Kendall's tau
           },
       {stop("wrong method")})

}

## ==== Minimum distance estimation ============================================

##' Distances for minimum distance estimation
##'
##' @title Distances for minimum distance estimation
##' @param u matrix of realizations (ideally) following U[0,1]^(d-1) or U[0,1]^d 
##' @param method distance methods available:
##'        mde.normal.CvM  = map to a chi-square distribution (Cramer-von Mises distance)
##'        mde.normal.KS   = map to a chi-square distribution (Kolmogorov-Smirnov distance)
##'        mde.log.CvM     = map to an Erlang distribution (Cramer-von Mises distance)
##'        mde.log.KS      = map to an Erlang distribution (Kolmogorov-Smirnov distance)
##' @return distance 
##' @author Marius Hofert
emde.dist <- function(u, method = c("mde.normal.CvM", "mde.normal.KS", "mde.log.CvM",
                         "mde.log.KS")){
    if(!is.matrix(u)) u <- rbind(u)
    d <- ncol(u)    
    n <- nrow(u)
    method <- match.arg(method) # match argument method
    switch(method,
           mde.normal.CvM = { # map to a chi-square distribution
               y <- sort(rowSums(qnorm(u)^2))
               Fvals <- pchisq(y, d)
               weights <- (2*(1:n)-1)/(2*n)
               1/(12*n) + sum((weights - Fvals)^2)
           },
           mde.normal.KS = { # map to a chi-square distribution
               y <- sort(rowSums(qnorm(u)^2))
               Fvals <- pchisq(y, d)
               i <- 1:n
               max(Fvals[i]-(i-1)/n, i/n-Fvals[i])
           },
           mde.log.CvM = { # map to an Erlang distribution
               y <- sort(rowSums(-log(u)))
               Fvals <- pgamma(y, shape = d) 
               weights <- (2*(1:n)-1)/(2*n)
               1/(12*n) + sum((weights - Fvals)^2)
           },
           mde.log.KS = { # map to an Erlang distribution
               y <- rowSums(-log(u))
               Fvals <- pgamma(y, shape = d) 
               i <- 1:n
               max(Fvals[i]-(i-1)/n, i/n-Fvals[i])
           },
           ## Note: The following multivariate distances turned out to be (far) too slow
           ## mde.SB = { # S_n^{(B)} from Genest et al. (2009)
           ##     sum1 <- sum(apply(1-u^2,1,prod))/2^(d-1)
           ##     f <- function(i,j) prod(1-apply(rbind(u[i,],u[j,]), 2, max)) 
           ##     sum2 <- 0
           ##     for(i in 1:n) sum2 <- sum2 + sum(unlist(lapply(1:n, function(j) f(i,j))))
           ##     n/3^d-sum1+sum2/n
           ## },
           ## mde.SC = { # S_n^{(C)} from Genest et al. (2009)
           ##     C.hat <- function(u.vec,u.mat) mean(apply(t(apply(u.mat, 1, 
           ##                                                       function(x) x <= 
           ##                                                       u.vec)), 1, all))
           ##     sum((apply(u, 1, C.hat, u.mat = u) - apply(u, 1, prod))^2)
           ## },
           stop("wrong distance method"))
}

##' Minimum distance estimation for nested Archimedean copulas
##'
##' @title Minimum distance estimation for nested Archimedean copulas
##' @param u matrix of realizations following the copula
##' @param cop outer_nacopula to be estimated
##' @param interval bivariate vector denoting the interval where optimization takes
##'        place
##' @param include.K boolean indicating whether the last component, K, is also used or not
##' @param method distance methods available, see emde.dist
##' @param ... additional parameters for optimize
##' @return minimum distance estimator; return value of optimize
##' @author Marius Hofert
emde <- function(u, cop, method = c("mde.normal.CvM", "mde.normal.KS", "mde.log.CvM",
                         "mde.log.KS"),
                 interval = paraOptInterval(u, cop@copula@name), 
                 include.K = ncol(u)<=5, ...)
{
    stopifnot(is(cop, "outer_nacopula"), is.numeric(d <- ncol(u)), d >= 1)
    if(length(cop@childCops))
        stop("currently, only Archimedean copulas are provided")
    method <- match.arg(method) # match argument method
    distance <- function(theta){ # distance to be minimized
        cop@copula@theta <- theta
        u. <- gnacopulatrafo(u, cop, n.MC=0, do.pseudo=FALSE, include.K=include.K) # transform data [don't use MC here; too slow]
        emde.dist(u., method)
    }
    optimize(distance, interval = interval, ...)	
}

## ==== Diagonal maximum likelihood estimation =================================

##' Density of the diagonal of a nested Archimedean copula
##'
##' @title Diagonal density of a nested Archimedean copula
##' @param u evaluation point
##' @param cop outer_nacopula to be estimated
##' @param log if TRUE the log-density is evaluated
##' @return density of the diagonal of cop
##' @author Marius Hofert
dDiag <- function(u, cop, log = FALSE){
    stopifnot(is(cop, "outer_nacopula"), all(0 <= u, u <= 1), (d <- max(cop@comp)) > 0)
    if(length(cop@childCops))
        stop("currently, only Archimedean copulas are provided")
    acop <- cop@copula
    th <- acop@theta
    if(log){
        log(d) + acop@psiDabs(d*acop@psiInv(u,th), th, log = TRUE) +
            acop@psiInvD1abs(u, th, log = TRUE)
    }else{
        d*acop@psiDabs(d*acop@psiInv(u,th), th)*(acop@psiInvD1abs(u,th))
    }
}

##' Maximum likelihood estimation based on the diagonal of a nested Archimedean copula
##'
##' @title Maximum likelihood estimation based on the diagonal of a nested Archimedean copula
##' @param u matrix of realizations following a copula
##' @param cop outer_nacopula to be estimated
##' @param interval bivariate vector denoting the interval where optimization takes
##'        place
##' @param ... additional parameters for optimize
##' @return diagonal maximum likelihood estimator; return value of optimize
##' @author Marius Hofert
edmle <- function(u, cop, interval = paraOptInterval(u, cop@copula@name), ...)
{
    stopifnot(is(cop, "outer_nacopula"), is.numeric(d <- ncol(u)), d >= 1, 
              max(cop@comp) == d) # dimension
    if(length(cop@childCops))
        stop("currently, only Archimedean copulas are provided")
    x <- apply(u,1,max) # data from the diagonal
    ## explicit estimator for Gumbel
    if(cop@copula@name == "Gumbel"){
	list(minimum = log(d)/(log(length(x))-log(sum(-log(x)))), objective = 0) # return value of the same structure as for optimize
    }else{
	## FIXME: clean-up unused code (or use method-switch for optim(x))
	## tau.hat.G <- copGumbel@tau(dmle.hat.G)
	##         start <-
	##             if(cop@name == "AMH" && tau.hat.G >= 0.33) {
	##                 ## AMH cannot attain a tau >= 1/3
	##                 0.99 # most likely, this does not converge then (especially if tau.hat.G >> 0.33)
	##             } else {
	##                 cop@tauInv(tau.hat.G) # convert tau.hat.G to the parameter of cop
	##             }
        ## optimize
	mLogL <- function(theta){ # -log-likelihood
            cop@copula@theta <- theta
            -sum(dDiag(x, cop=cop, log=TRUE)) 
        }
	optimize(mLogL, interval=interval, ...)
	##optimx::optimx(start, function(theta) # -log-Likelihood of the diagonal
        ##               sum(cop@dDiag(x,theta,log = TRUE)),
        ##               lower= min(cop@paraSubInterval),
        ##               upper= max(cop@paraSubInterval),
        ##               ## For box constraints ('lower', 'upper'), method must be one of
        ##               ##  c("L-BFGS-B", "nlminb", "spg", "Rcgmin", "Rvmmin", "bobyqa")
        ##               method = "bobyqa",
        ##              ...)
    }
}

## ==== (Simulated) maximum likelihood estimation ==============================

##' (Simulated) maximum likelihood estimation for nested Archimedean copulas
##'
##' @title (Simulated) maximum likelihood estimation for nested Archimedean copulas
##' @param u matrix of realizations following the copula
##' @param cop outer_nacopula to be estimated
##' @param n.MC if > 0 SMLE is applied with sample size equal to n.MC; otherwise,
##'        MLE is applied
##' @param interval bivariate vector denoting the interval where optimization takes
##'        place (with default given by the slot paraSubInterval)
##' @param ... additional parameters for optimize
##' @return (simulated) maximum likelihood estimator; return value of optimize
##' @author Marius Hofert
emle <- function(u, cop, n.MC=0, interval=paraOptInterval(u, cop@copula@name), ...)
{
    stopifnot(is(cop, "outer_nacopula"))
    if(length(cop@childCops))
	stop("currently, only Archimedean copulas are provided")
    ## compute initial value based on either etau or edmle
    ## FIXME: clean-up unused code
    ## this function used to have the argument "initial = c("tau.mean", "theta.mean","diag")"
    ## with the documentation
    ##     initial method used for finding an initial value (either etau with
    ## 	   option "tau.mean" or "theta.mean" (denoted by "tau.mean" and "tau.mean",
    ##	   respectively), or edmle (denoted by "diag"))
    ## and the following code (here)
    ## start <- switch(match.arg(initial),
    ##                     tau.mean   = etau (u,acop, method="tau.mean"),
    ##                     theta.mean = etau (u,acop, method="theta.mean"),
    ##                     diag       = edmle(u,acop,...),
    ##                     stop("wrong argument initial"))
    ## optimize
    mLogL <- function(theta){ # -log-likelihood
        cop@copula@theta <- theta
        -sum(dnacopula(cop, u, n.MC=n.MC, log=TRUE))
    }
    optimize(mLogL, interval=interval, ...)
    ## optimx::optimx(start, mLogL,
    ##                    lower= min(acop@paraSubInterval),
    ##                    upper= max(acop@paraSubInterval),
    ##                    ## For box constraints ('lower', 'upper'), method must be one of
    ##                    ##  c("L-BFGS-B", "nlminb", "spg", "Rcgmin", "Rvmmin", "bobyqa")
    ##                    method = "bobyqa",
    ##                    ...)

}

## ==== Estimation wrapper =====================================================

##' Computes the pseudo-observations for the given data matrix
##'
##' @title Pseudo-observations 
##' @param x matrix of random variates to be converted to pseudo-observations
##' @return pseudo-observations (matrix of the same dimensions as x)
##' @author Marius Hofert
pobs <- function(x) apply(x,2,rank)/(nrow(x)+1)

##' Computes different parameter estimates for a nested Archimedean copula
##'
##' @title Estimation procedures for nested Archimedean copulas
##' @param x data matrix
##' @param cop outer_nacopula to be estimated
##' @param method estimation method; can be
##'        "mle"             MLE
##'        "smle"            SMLE
##'        "dmle"            MLE based on the diagonal
##'        "mde.normal.CvM"  minimum distance estimation based on the chisq distribution and CvM distance
##'        "mde.normal.KS"   minimum distance estimation based on the chisq distribution and KS distance
##'        "mde.log.CvM"     minimum distance estimation based on the Erlang distribution and CvM distance
##'        "mde.log.KS"      minimum distance estimation based on the Erlang distribution and KS distance
##'        "tau.tau.mean"    averaged pairwise Kendall's tau estimator
##'        "tau.theta.mean"  average of Kendall's tau estimators
##'        "beta"            multivariate Blomqvist's beta estimator
##' @param n.MC if > 0 it denotes the sample size for SMLE
##' @param interval initial optimization interval for "mle", "smle", and "dmle" (i.e., emle, dmle)
##' @param do.pseudo  logical indicating if 'x' should be "pseudo-transformed"
##' @param xargs additional arguments for the estimation procedures
##' @param ... additional parameters for optimize
##' @return estimated value/vector according to the chosen method
##' @author Marius Hofert
enacopula <- function(x, cop, method=c("mle", "smle", "dmle", "mde.normal.CvM", 
                              "mde.normal.KS", "mde.log.CvM", "mde.log.KS",
                              "tau.tau.mean", "tau.theta.mean", "beta"), 
                      n.MC = if(method=="smle") 10000 else 0, 
                      interval=paraOptInterval(u, cop@copula@name), 
                      do.pseudo=FALSE, xargs=list(), ...)
{
    
    ## setup 
    stopifnot(is(cop, "outer_nacopula"), is.list(xargs))
    if(length(cop@childCops))
        stop("currently, only Archimedean copulas are provided")
    u <- if(do.pseudo) pobs(x) else x
    method <- match.arg(method) 

    ## main part
    res <- switch(method,
                  mle =            do.call(emle, c(list(u, cop, 
                  interval = interval, ...), xargs)),
                  smle =           do.call(emle, c(list(u, cop, n.MC = n.MC, 
                  interval = interval, ...), xargs)),
                  dmle =           do.call(edmle, c(list(u, cop, 
                  interval = interval, ...), xargs)),
                  mde.normal.CvM = do.call(emde, c(list(u, cop, "mde.normal.CvM", 
                  interval = interval, ...), xargs)),
                  mde.normal.KS =  do.call(emde, c(list(u, cop, "mde.normal.KS", 
                  interval = interval, ...), xargs)),
                  mde.log.CvM =    do.call(emde, c(list(u, cop, "mde.log.CvM", 
                  interval = interval, ...), xargs)), 
                  mde.log.KS =     do.call(emde, c(list(u, cop, "mde.log.KS", 
                  interval = interval, ...), xargs)),
                  tau.tau.mean =   do.call(etau, c(list(u, cop, "tau.mean", ...), 
                  xargs)),
                  tau.theta.mean = do.call(etau, c(list(u, cop, "theta.mean", ...), 
                  xargs)),
                  beta =           do.call(ebeta, c(list(u, cop, 
                  interval = interval, ...), xargs)),
                  stop("wrong estimation method"))

    ## FIXME: deal with result, check details, give warnings

    ## return the estimate
    switch(method,
           mle =            res$minimum,
           smle =           res$minimum,
           dmle =           res$minimum,
           mde.normal.CvM = res$minimum,
           mde.normal.KS =  res$minimum, 
           mde.log.CvM =    res$minimum,  
           mde.log.KS =     res$minimum,
           tau.tau.mean =   res,
           tau.theta.mean = res,
           beta =           res$root,
           stop("wrong estimation method"))

}
