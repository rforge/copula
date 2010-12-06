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
##' @param u data
##' @param family Archimedean family
##' @param h for enlarging the tau-interval
##' @return initial interval which can be used for optimization (e.g., for emle)
##' @author Marius Hofert
paraOptInterval <- function(u, family, h = 0.15){
    theta.hat.G <- edmle(u, getAcop("Gumbel"))$minimum
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

##' Compute the sample version of Blomqvist's beta for an Archimedean copula,
##' see, e.g., Schmid and Schmidt (2007) "Nonparametric inference on multivariate
##' versions of Blomqvist’s beta and related measures of tail dependence"
##' @param u matrix of realizations following the copula
##' @param scaling if FALSE then the scaling factors 2^(d-1)/(2^(d-1)-1) and
##'                2^(1-d) are omitted
##' @return sample version of multivariate Blomqvist beta
##' @author Marius Hofert
beta.hat <- function(u, scaling = TRUE){
    stopifnot(all(0 <= u, u < 1))
    less.u <- u <= 0.5
    prod1 <- apply( less.u, 1, all)
    prod2 <- apply(!less.u, 1, all)
    sum.prod <- prod1 + prod2
    b <- mean(sum.prod)
    d <- ncol(u)
    if(scaling) {T <- 2^(d-1); (T*b - 1)/(T - 1)} else b
}

##' Compute the population version of Blomqvist's beta for an Archimedean copula
##' @param cop acopula to be estimated
##' @param theta copula parameter
##' @param d dimension
##' @param scaling if FALSE then the scaling factors 2^(d-1)/(2^(d-1)-1) and
##'                2^(1-d) are omitted
##' @return population version of multivariate Blomqvist beta
##' @author Marius Hofert & Martin Maechler
beta. <- function(cop, theta, d, scaling = TRUE) {
    stopifnot(cop@paraConstr(theta)) # also checks length(theta)
    j <- 1:d
    b <- pnacopula(onacopulaL(cop@name, list(theta, j)), u = rep.int(0.5, d)) +1
    if(d < 30) { ## for small d, using log(.) and exp(.) is wasteful (and looses precision)
	diag. <- function(k)
	    pnacopula(onacopulaL(cop@name, list(theta,seq_len(k))),
		      u = rep.int(0.5, k))
	diags <- unlist(lapply(j, diag.))
	ch.dia <- choose(d,j) * diags

    } else { ## "large" d

	log.diag <- function(k)
	    log(pnacopula(onacopulaL(cop@name, list(theta,seq_len(k))),
			  u = rep.int(0.5, k)))
	ldiags <- unlist(lapply(j, log.diag))
	ch.dia <- exp(lchoose(d,j)+ldiags)
    }
    b <- b + sum((-1)^j * ch.dia)
    if(scaling) { T <- 2^(d-1); (T*b - 1)/(T - 1)} else b
}

##' Method-of-moment-like estimation of Archimedean copulas based on a
##' multivariate version of Blomqvist's beta
##' @param u matrix of realizations following the copula
##' @param cop acopula to be estimated
##' @param interval bivariate vector denoting the interval where optimization takes
##'        place
##' @param ... additional parameters for safeUroot
##' @return Blomqvist beta estimator; return value of safeUroot (more or less
##'	    equal to the return value of uniroot)
##' @author Marius Hofert
ebeta <- function(u,cop,interval = paraOptInterval(u, cop@name),...){
    ## Note: We do not need the constants 2^(d-1)/(2^(d-1)-1) and 2^(1-d) here,
    ##	     since we equate the population and sample versions of Blomqvist's
    ##       beta anyway.
    b.hat <- beta.hat(u, scaling = FALSE)
    d <- ncol(u)
    safeUroot(function(theta){beta.(cop,theta,d,scaling = FALSE) - b.hat},
              interval = interval, Sig = +1, check.conv = TRUE, ...)
}

## ==== Kendall's tau ==========================================================

##' Compute pairwise Kendall's tau estimator for Archimedean copulas
##' @param u matrix of realizations following the copula
##' @param cop acopula to be estimated
##' @param method tau.mean indicates that the average of the sample versions of
##'               Kendall's tau are computed first and then theta is determined;
##'               theta.mean stands for first computing all Kendall's tau
##'               estimators and then returning the mean of these estimators
##' @param ... additional arguments to cor()
##' @return averaged pairwise cor() estimators
##' @author Marius Hofert
etau <- function(u,cop,method = c("tau.mean","theta.mean"),...)
{
    stopifnot(is(cop,"acopula"))
    tau.hat.mat <- cor(u, method="kendall",...) # matrix of pairwise tau()
    tau.hat <- tau.hat.mat[upper.tri(tau.hat.mat)] # all tau hat's

    method <- match.arg(method)
    switch(method,
           "tau.mean" = {
               cop@tauInv(mean(tau.hat)) # Kendall's tau corresponding to the mean of the tau hat's
           },
           "theta.mean" = {
               mean(cop@tauInv(tau.hat)) # mean of the Kendall's tau
           },
       {stop("wrong method")})

}

## ==== Diagonal maximum likelihood estimation =================================

##' Density of the diagonal of an Archimedean copula
##' @param u evaluation point
##' @param cop acopula
##' @param d dimension
##' @param log if TRUE the log-density is evaluated
##' @return density of the diagonal of cop
##' @author Marius Hofert
dDiag <- function(u, cop, d, log = FALSE){
    stopifnot(is(cop,"acopula")) # dimension
    stopifnot(all(0 <= u, u <= 1))
    th <- cop@theta
    if(log){
        log(d) + cop@psiDabs(d*cop@psiInv(u,th), th, log = TRUE) +
            cop@psiInvD1abs(u, th, log = TRUE)
    }else{
        d*cop@psiDabs(d*cop@psiInv(u,th), th)*(cop@psiInvD1abs(u,th))
    }
}

##' Maximum likelihood estimation based on the diagonal of the Archimedean copula
##' @param u matrix of realizations following the copula
##' @param cop acopula to be estimated
##' @param interval bivariate vector denoting the interval where optimization takes
##'        place
##' @param ... additional parameters for optimize
##' @return diagonal maximum likelihood estimator; return value of optimize
##' @author Marius Hofert
edmle <- function(u, cop, interval = paraOptInterval(u, cop@name), ...)
{
    stopifnot(is(cop,"acopula"), is.numeric(d <- ncol(u)), d >= 1) # dimension
    x <- apply(u,1,max) # data from the diagonal
    ## explicit estimator for Gumbel
    if(cop@name == "Gumbel"){
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
            cop@theta <- theta
            -sum(dDiag(x,cop,d,log = TRUE)) ## FIXME: maybe only sum over those values that are finite and let the user know how many these are (?)
        }
	optimize(mLogL, interval = interval, ...)
	##optimx::optimx(start, function(theta) # -log-Likelihood of the diagonal
        ##               sum(cop@dDiag(x,theta,d,log = TRUE)),
        ##               lower= min(cop@paraSubInterval),
        ##               upper= max(cop@paraSubInterval),
        ##               ## For box constraints ('lower', 'upper'), method must be one of
        ##               ##  c("L-BFGS-B", "nlminb", "spg", "Rcgmin", "Rvmmin", "bobyqa")
        ##               method = "bobyqa",
        ##              ...)
    }
}

## ==== (Simulated) maximum likelihood estimation ==============================

##' (Simulated) maximum likelihood estimation for Nested Archimedean copulas
##' @param u matrix of realizations following the copula
##' @param cop nacopula to be estimated
##' @param MC if provided SMLE is applied with sample size equal to MC; otherwise,
##'        MLE is applied
##' @param interval bivariate vector denoting the interval where optimization takes
##'        place (with default given by the slot paraSubInterval)
##' @param ... additional parameters for optimize
##' @return (simulated) maximum likelihood estimator; return value of optimize
##' @author Marius Hofert
emle <- function(u, cop, MC, interval = paraOptInterval(u, cop@copula@name), ...)
{
    stopifnot(is(cop,"nacopula"))
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
    if(missing(MC)) MC <- NULL # set it to NULL to make the call from dnacopula easier
    ## optimize
    mLogL <- function(theta){ # -log-likelihood
        cop@copula@theta <- theta
### FIXME: maybe only sum over those values that are finite and let the user know how many these are (?)
        -sum(dnacopula(cop, u, MC = MC, log = TRUE))
    }
    optimize(mLogL, interval = interval, ...)
    ## optimx::optimx(start, mLogL,
    ##                    lower= min(acop@paraSubInterval),
    ##                    upper= max(acop@paraSubInterval),
    ##                    ## For box constraints ('lower', 'upper'), method must be one of
    ##                    ##  c("L-BFGS-B", "nlminb", "spg", "Rcgmin", "Rvmmin", "bobyqa")
    ##                    method = "bobyqa",
    ##                    ...)

}

## ==== Estimation wrapper =====================================================

##' Computes the (scaled) pseudo-observations for the given data matrix
##' @param x matrix of random variates to be converted to pseudo-observations
##' @return pseudo-observations (matrix of the same dimensions as x)
##' @author Marius Hofert
pobs <- function(x) apply(x,2,rank)/(nrow(x)+1)

##' Computes different parameter estimates for a nested Archimedean copula
##' @param x data matrix
##' @param cop nacopula to be estimated
##' @param method estimation method; can be
##'        "mle"             MLE
##'        "smle"            SMLE
##'        "tau.tau.mean"    averaged pairwise Kendall's tau estimator
##'        "tau.theta.mean"  average of Kendall's tau estimators
##'        "dmle"            MLE based on the diagonal
##'        "beta"            multivariate Blomqvist's beta estimator
##' FIXME: the docu (and enacopula) used to contain:
##'        "mle.tau.mean"     MLE with initial value found by averaged pairwise Kendall's tau
##'        "mle.theta.mean"   MLE with initial value found by Kendall's tau estimators averaged
##'        "mle.diag"         MLE with initial value found via DMLE
##'        "smle.tau.mean"    SMLE with initial value found by averaged pairwise Kendall's tau
##'        "smle.theta.mean"  SMLE with initial value found by Kendall's tau estimators averaged
##'        "smle.diag"        SMLE with initial value found via DMLE
##' @param MC if provided (and not NULL) MC denotes the sample size for SMLE
##' @param interval initial optimization interval for "mle", "smle", and "dmle" (i.e., emle, dmle)
##' @param do.pseudo  logical indicating if 'x' should be "pseudo-transformed"
##' @param ... additional parameters for optimize
##' @return estimated value/vector according to the chosen method
##' @author Marius Hofert
enacopula <- function(x, cop, method = c("mle","smle","tau.tau.mean","tau.theta.mean",
                              "dmle","beta"), MC, interval =
                      paraOptInterval(u, cop@copula@name), do.pseudo = FALSE, ...)
{
    ## setup cop
    stopifnot(is(cop, "outer_nacopula"))
    if(length(cop@childCops))
        stop("currently, only Archimedean copulas are provided")
    acop <- cop@copula

    ## setup x
    u <- if(do.pseudo) pobs(x) else x

    method <- match.arg(method)

    ## check if MC is given for SMLE
    if(mMC <- missing(MC)) MC <- NULL
    if(mMC && method == "smle") stop("smle needs the sample size MC")

    ## main part
    res <- switch(method,
                  mle =            emle (u,  cop, interval = interval,...),
                  smle =           emle (u,  cop, MC = MC, interval = interval,...),
                  tau.tau.mean =   etau (u, acop, "tau.mean", ...),
                  tau.theta.mean = etau (u, acop, "theta.mean", ...),
                  dmle =           edmle(u, acop, interval = interval,...),
                  beta =           ebeta(u, acop, interval = interval,...),
                  stop("wrong estimation method"))

    ## FIXME: deal with result, check details, give warnings

    ## return the estimate
    switch(method,
           mle =            res$minimum,
           smle =           res$minimum,
           tau.tau.mean =   res,
           tau.theta.mean = res,
           dmle =           res$minimum,
           beta =           res$root,
           stop("wrong estimation method"))

}
