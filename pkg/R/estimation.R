## Copyright (C) 2010--2011 Marius Hofert and Martin Maechler
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

## ==== initial interval/value for optimization procedures =====================

##' Compute an initial interval or value for optimization/estimation routines
##' (only a heuristic; if this fails, choose your own interval or value)
##'
##' @title Compute an initial interval or value for estimation procedures
##' @param family Archimedean family
##' @param tau.range vector containing lower and upper admissible Kendall's tau
##' @param value logical determining if an interval (the default) or an initial
##'        value should be returned
##' @param u matrix of realizations following a copula
##' @param method method for obtaining initial values
##' @param ... further arguments to cor() for method="tau.mean"
##' @return initial interval or value which can be used for optimization
##' @author Marius Hofert
initOpt <- function(family, tau.range=NULL, value=FALSE, u, method=c("tau.Gumbel",
                                                            "tau.mean"), ...){
    cop <- getAcop(family)
    if(is.null(tau.range)){
	eps <- 1e-8
        tau.range <- switch(cop@name, # limiting (attainable) taus that can be dealt with in estimation/optimization/root-finding
                            "AMH" = { c(0, 1/3-eps) },
                            "Clayton" = { c(eps, 0.95) },
                            "Frank" = { c(eps, 0.94) }, # beyond that, estimation.gof() fails for ebeta()!
                            "Gumbel" = { c(0, 0.95) },
                            "Joe" = { c(0, 0.95) },
                            stop("unsupported family for initOpt"))
    }
    if(!value) return(cop@tauInv(tau.range)) # for value=FALSE, u is not required
    stopifnot(length(dim(u)) == 2)
    method <- match.arg(method)
    ## estimate Kendall's tau
    tau.hat <- switch(method,
                      "tau.Gumbel"={
                          x <- apply(u, 1, max)
                          theta.hat.G <- log(ncol(u))/(log(length(x))-log(sum(-log(x)))) # direct formula from edmle for Gumbel
                          copGumbel@tau(theta.hat.G)
                      },
                      "tau.mean"={
                          tau.hat.mat <- cor(u, method="kendall", ...) # matrix of pairwise tau()
                          mean(tau.hat.mat[upper.tri(tau.hat.mat)]) # mean of estimated taus
                      },
                      stop("wrong method for initOpt"))
    ## truncate if required
    if(tau.hat < tau.range[1]){
        cop@tauInv(tau.range[1])
    }else if(tau.hat > tau.range[2]){
        cop@tauInv(tau.range[2])
    }else{
        cop@tauInv(tau.hat)
    }
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
beta.hat <- function(u, scaling = FALSE) {
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
beta. <- function(cop, theta, d, scaling=FALSE) {
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
ebeta <- function(u, cop, interval=initOpt(cop@copula@name), ...) {
    stopifnot(is(cop, "outer_nacopula"), is.numeric(d <- ncol(u)), d >= 2,
              max(cop@comp) == d)
    if(length(cop@childCops))
        stop("currently, only Archimedean copulas are provided")
    ## Note: We do not need the constants 2^(d-1)/(2^(d-1)-1) and 2^(1-d) here,
    ##	     since we equate the population and sample versions of Blomqvist's
    ##       beta anyway.
    b.hat <- beta.hat(u, scaling = TRUE)
    d <- ncol(u)
    safeUroot(function(theta) {beta.(cop@copula, theta, d, scaling=TRUE) - b.hat},
              interval=interval, Sig=+1, check.conv=TRUE, ...)
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
##' @param warn logical determining if warnings are produced (for AMH and in
##'        general for pairwise sample versions of Kendall's tau < 0)
##' @param ... additional arguments to cor()
##' @return averaged pairwise cor() estimators
##' @author Marius Hofert
etau <- function(u, cop, method = c("tau.mean", "theta.mean"), warn=TRUE, ...)
{
    stopifnot(is(cop, "outer_nacopula"), is.numeric(d <- ncol(u)), d >= 2,
              max(cop@comp) == d)
    if(length(cop@childCops))
        stop("currently, only Archimedean copulas are provided")
    tau.hat.mat <- cor(u, method="kendall",...) # matrix of pairwise tau()
    tau.hat <- tau.hat.mat[upper.tri(tau.hat.mat)] # all tau hat's
    ## check if there are tau's < 0
    if((l <- length(ind <- which(tau.hat < 0))) > 0){
	if(warn){
	    ct <- sort(tau.hat[ind])
            ct <- if(l > 3) paste(paste(format(ct[1:3]), collapse=", "),
                                  ", ...",sep="") else format(ct)
            warning("etau: There ",if(l==1) "is " else "are ",l," pairwise Kendall's tau < 0, ",
                    if(l==1) "" else "the smallest ","being tau = ", ct)
        }
        tau.hat[ind] <- 0
    }
    ## apply tauInv in the appropriate way
    tau_inv <- if(cop@copula@name == "AMH")
	function(tau) cop@copula@tauInv(tau, check=FALSE, warn=warn) else cop@copula@tauInv
    method <- match.arg(method)
    switch(method,
           "tau.mean" = {
               tau_inv(mean(tau.hat)) # Kendall's tau corresponding to the mean of the tau hat's
           },
           "theta.mean" = {
               mean(tau_inv(tau.hat)) # mean of the Kendall's tau
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
                         "mde.log.KS")) {
    if(!is.matrix(u)) u <- rbind(u)
    d <- ncol(u)
    n <- nrow(u)
    method <- match.arg(method) # match argument method
    switch(method,
           "mde.normal.CvM" = { # map to a chi-square distribution
               y <- sort(rowSums(qnorm(u)^2))
               Fvals <- pchisq(y, d)
               weights <- (2*(1:n)-1)/(2*n)
               1/(12*n) + sum((weights - Fvals)^2)
           },
           "mde.normal.KS" = { # map to a chi-square distribution
               y <- sort(rowSums(qnorm(u)^2))
               Fvals <- pchisq(y, d)
               i <- 1:n
               max(Fvals[i]-(i-1)/n, i/n-Fvals[i])
           },
           "mde.log.CvM" = { # map to an Erlang distribution
               y <- sort(rowSums(-log(u)))
               Fvals <- pgamma(y, shape = d)
               weights <- (2*(1:n)-1)/(2*n)
               1/(12*n) + sum((weights - Fvals)^2)
           },
           "mde.log.KS" = { # map to an Erlang distribution
               y <- rowSums(-log(u))
               Fvals <- pgamma(y, shape = d)
               i <- 1:n
               max(Fvals[i]-(i-1)/n, i/n-Fvals[i])
           },
           ## Note: The following multivariate distances turned out to be (far) too slow
           ## "mde.SB" = { # S_n^{(B)} from Genest et al. (2009)
           ##     sum1 <- sum(apply(1-u^2,1,prod))/2^(d-1)
           ##     f <- function(i,j) prod(1-apply(rbind(u[i,],u[j,]), 2, max))
           ##     sum2 <- 0
           ##     for(i in 1:n) sum2 <- sum2 + sum(unlist(lapply(1:n, function(j) f(i,j))))
           ##     n/3^d-sum1+sum2/n
           ## },
           ## "mde.SC" = { # S_n^{(C)} from Genest et al. (2009)
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
                 interval = initOpt(cop@copula@name),
                 include.K = ncol(u)<=5, ...)
{
    stopifnot(is(cop, "outer_nacopula"), is.numeric(d <- ncol(u)), d >= 2,
              max(cop@comp) == d)
    if(length(cop@childCops))
        stop("currently, only Archimedean copulas are provided")
    method <- match.arg(method) # match argument method
    distance <- function(theta) { # distance to be minimized
        cop@copula@theta <- theta
        u. <- gnacopulatrafo(u, cop, n.MC=0, include.K=include.K) # transform data [don't use MC here; too slow]
        emde.dist(u., method)
    }
    optimize(distance, interval = interval, ...)
}

## ==== Diagonal maximum likelihood estimation =================================

##' Density of the diagonal of a nested Archimedean copula
##'
##' @title Diagonal density of a nested Archimedean copula
##' @param u evaluation point in [0,1]
##' @param cop outer_nacopula
##' @param log if TRUE the log-density is evaluated
##' @return density of the diagonal of cop
##' @author Marius Hofert
dDiag <- function(u, cop, log=FALSE) {
    stopifnot(is(cop, "outer_nacopula"), (d <- max(cop@comp)) >= 2)
    if(length(cop@childCops)) {
        stop("currently, only Archimedean copulas are provided")
    }
    else cop@copula@dDiag(u, theta=cop@copula@theta, d=d, log=log)
}

##' @title Generic density of the diagonal of d-dim. Archimedean copula
##' @param u evaluation point in [0, 1]
##' @param d dimension
##' @param cop acopula
##' @param log if TRUE the log-density is evaluated
##' @return density of the diagonal of cop
##' @author Martin Maechler
dDiagA <- function(u, d, cop, log=FALSE) {
    stopifnot(is.finite(th <- cop@theta), d >= 2)
    if(log) {
        log(d) + cop@psiDabs(d*cop@psiInv(u, th), th, log=TRUE) +
            cop@psiInvD1abs(u, th, log=TRUE)
    } else {
        d * cop@psiDabs(d*cop@psiInv(u, th), th) * cop@psiInvD1abs(u, th)
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
edmle <- function(u, cop, interval=initOpt(cop@copula@name), ...)
{
    stopifnot(is(cop, "outer_nacopula"), is.numeric(d <- ncol(u)), d >= 2,
              max(cop@comp) == d) # dimension
    if(length(cop@childCops))
        stop("currently, only Archimedean copulas are provided")
    x <- apply(u, 1, max) # data from the diagonal
    ## explicit estimator for Gumbel
    if(cop@copula@name == "Gumbel") {
	list(minimum = log(d)/(log(length(x))-log(sum(-log(x)))), objective = 0) # return value of the same structure as for optimize
    } else {
        ## optimize
	nlogL <- function(theta) # -log-likelihood of the diagonal
            -sum(cop@copula@dDiag(x, theta=theta, d=d, log=TRUE))
	optimize(nlogL, interval=interval, ...)
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
##'        place
##' @param ... additional parameters for optimize
##' @return (simulated) maximum likelihood estimator; return value of optimize
##' @author Marius Hofert
##' note: this is a *fast* version (based on optimize()) called from enacopula
.emle <- function(u, cop, n.MC=0, interval=initOpt(cop@copula@name), ...)
{
    stopifnot(is(cop, "outer_nacopula"))
    if(length(cop@childCops))
	stop("currently, only Archimedean copulas are provided")
    ## optimize
    mLogL <- function(theta) { # -log-likelihood
        cop@copula@theta <- theta
        -sum(dnacopula(cop, u, n.MC=n.MC, log=TRUE))
    }
    optimize(mLogL, interval=interval, ...)
}


##' (Simulated) maximum likelihood estimation for nested Archimedean copulas
##'
##' @title (Simulated) maximum likelihood estimation for nested Archimedean copulas
##' @param u matrix of realizations following the copula
##' @param cop outer_nacopula to be estimated
##' @param n.MC if > 0 SMLE is applied with sample size equal to n.MC; otherwise,
##'        MLE is applied
##' @param interval bivariate vector denoting the interval where optimization takes
##'        place
##' @param ... additional parameters for optimize
##' @return an "mle2" object with the (simulated) maximum likelihood estimator.
##' @author Martin Maechler and Marius Hofert
##' note: this is the *slower* version which also allows for profiling
emle <- function(u, cop, n.MC=0, optimizer="optimize", method,
		 interval=initOpt(cop@copula@name),
                 ##vvv awkward to be needed, but it is - by mle2():
                 start = list(theta=initOpt(cop@copula@name, value=TRUE, u=u)),
                 ...)
{
    stopifnot(is(cop, "outer_nacopula"), is.numeric(d <- ncol(u)), d >= 2,
              max(cop@comp) == d)
    ## nLL <- function(theta) { # -log-likelihood
    ##	   cop@copula@theta <- theta
    ##	   -sum(dnacopula(cop, u, n.MC=n.MC, log=TRUE))
    ## }
    if(length(cop@childCops))
	stop("currently, only Archimedean copulas are provided")
    else ## For (*non*-nested) copulas only:
	nLL <- function(theta)  # -(log-likelihood)
	    -sum(cop@copula@dacopula(u, theta, n.MC=n.MC, log=TRUE))

    ## optimization
    if(!(is.null(optimizer) || is.na(optimizer))) {
        stopifnot(require("bbmle"))
	if(optimizer == "optimize")
	    mle2(minuslogl = nLL, optimizer = "optimize",
		 lower = interval[1], upper = interval[2],
		 ##vvv awkward to be needed, but it is - by mle2():
		 start=start, ...)
	else if(optimizer == "optim") {
	    message(" optimizer = \"optim\" -- using mle2(); consider optimizer=NULL instead")
	    mle2(minuslogl = nLL, optimizer = "optim", method = method,
		 start=start, ...)
	}
	else ## "general"
	    mle2(minuslogl = nLL, optimizer = optimizer,
		 ##vvv awkward to be needed, but it is - by mle2():
		 start=start, ...)
    }
    else
	## use optim() .. [which uses suboptimal method for 1D, but provides Hessian]
	mle(minuslogl = nLL, method = method, start=start, ...)
}

## ==== Estimation wrapper =====================================================

##' Computes the pseudo-observations for the given data matrix
##'
##' @title Pseudo-observations
##' @param x matrix of random variates to be converted to pseudo-observations
##' @return pseudo-observations (matrix of the same dimensions as x)
##' @author Marius Hofert
pobs <- function(x, na.last="keep", ties.method=c("average", "first", "random", "max", "min"))
    apply(x,2,rank, na.last=na.last, ties.method=ties.method)/(nrow(x)+1)

##' Computes different parameter estimates for a nested Archimedean copula
##'
##' @title Estimation procedures for nested Archimedean copulas
##' @param u data matrix
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
##' @param interval initial optimization interval for "mle", "smle", and "dmle"
##' @param xargs additional arguments for the estimation procedures
##' @param ... additional parameters for optimize
##' @return estimated value/vector according to the chosen method
##' @author Marius Hofert
enacopula <- function(u, cop, method=c("mle", "smle", "dmle", "mde.normal.CvM",
                              "mde.normal.KS", "mde.log.CvM", "mde.log.KS",
                              "tau.tau.mean", "tau.theta.mean", "beta"),
                      n.MC = if(method=="smle") 10000 else 0,
                      interval=initOpt(cop@copula@name),
                      xargs=list(), ...)
{

    ## setup
    if(!is.matrix(u)) u <- rbind(u)
    stopifnot(0 <= u, u <= 1, is(cop, "outer_nacopula"), (d <- ncol(u)) >= 2,
              max(cop@comp) == d, n.MC >= 0, is.list(xargs))
    if(length(cop@childCops))
        stop("currently, only Archimedean copulas are provided")
    method <- match.arg(method)

    ## main part
    res <- switch(method,
                  "mle" =            do.call(.emle, c(list(u, cop,
                  interval = interval, ...), xargs)),
                  "smle" =           do.call(.emle, c(list(u, cop, n.MC = n.MC,
                  interval = interval, ...), xargs)),
                  "dmle" =           do.call(edmle, c(list(u, cop,
                  interval = interval, ...), xargs)),
                  "mde.normal.CvM" = do.call(emde, c(list(u, cop, "mde.normal.CvM",
                  interval = interval, ...), xargs)),
                  "mde.normal.KS" =  do.call(emde, c(list(u, cop, "mde.normal.KS",
                  interval = interval, ...), xargs)),
                  "mde.log.CvM" =    do.call(emde, c(list(u, cop, "mde.log.CvM",
                  interval = interval, ...), xargs)),
                  "mde.log.KS" =     do.call(emde, c(list(u, cop, "mde.log.KS",
                  interval = interval, ...), xargs)),
                  "tau.tau.mean" =   do.call(etau, c(list(u, cop, "tau.mean", ...),
                  xargs)),
                  "tau.theta.mean" = do.call(etau, c(list(u, cop, "theta.mean", ...),
                  xargs)),
                  "beta" =           do.call(ebeta, c(list(u, cop,
                  interval = interval, ...), xargs)),
                  stop("wrong estimation method for enacopula"))

    ## FIXME: deal with result, check details, give warnings

    ## return the estimate
    switch(method,
           "mle" =            res$minimum,
           "smle" =           res$minimum,
           "dmle" =           res$minimum,
           "mde.normal.CvM" = res$minimum,
           "mde.normal.KS" =  res$minimum,
           "mde.log.CvM" =    res$minimum,
           "mde.log.KS" =     res$minimum,
           "tau.tau.mean" =   res,
           "tau.theta.mean" = res,
           "beta" =           res$root,
           stop("wrong estimation method"))

}
