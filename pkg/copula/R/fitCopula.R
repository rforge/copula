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


### Class and methods for fitCopula ############################################

setClass("fitCopula",
	 representation(copula = "copula"),
	 contains="fittedMV" #-> ./Classes.R
	 ## FIXME , validity = function(object) TRUE
	 )

setMethod("paramNames", "fitCopula", function(x) x@copula@param.names)

print.fitCopula <- function(x, digits = max(3, getOption("digits") - 3),
                            signif.stars = getOption("show.signif.stars"), ...)
{
    sfit <- summary.fitCopula(x)
    cat("fitCopula() estimation based on '", x@method, "'\nand a sample of size ",
	x@nsample, ".\n", sep="")
    printCoefmat(sfit$coefficients, digits = digits, signif.stars = signif.stars,
		 na.print = "NA", ...)
    if (!is.na(sfit$loglik))
	cat("The maximized loglikelihood is ", format(sfit$loglik, digits=digits), "\n")
    if (!is.na(sfit$convergence)) {
	if(sfit$convergence)
            cat("Convergence problems: code is", sfit$convergence, "see ?optim.\n")
	else cat("Optimization converged\n")
    }
    if(!is.null(cnts <- x@fitting.stats$counts) && !all(is.na(cnts))) {
	cat("Number of loglikelihood evaluations:\n"); print(cnts, ...)
    }
    ## FIXME:  Show parts of the @copula  slot !!!!!!
    invisible(x)
}

## FIXME:  print( summary( fitC ) )  gives *less*  than print ( fitC ) !!
summary.fitCopula <- function(object, ...) {
    estimate <- object@estimate
    se <- sqrt(diag(object@var.est))
    zvalue <- estimate / se
    pvalue <- 2* pnorm(abs(zvalue), lower.tail=FALSE)
    coef <- cbind(Estimate = estimate, "Std. Error" = se,
                  "z value" = zvalue, "Pr(>|z|)" = pvalue)
    rownames(coef) <- paramNames(object)
    structure(class = "summary.fitCopula",
              list(method = object@method,
                   loglik = object@loglik,
                   convergence = object@fitting.stats[["convergence"]],
                   coefficients = coef))
}

setMethod("show", signature("fitCopula"),
	  function(object) print.fitCopula(object))


### Auxiliary functions ########################################################

##' @title Construct Parameter Matrix from Matrix of Kendall's Taus
##' @param cop Copula object (typically with NA parameters)
##' @param x The data in [0,1]^d or IR^d
##' @param method The rank correlation measure used
##' @param posDef A logical indicating whether a proper correlation matrix
##'        is computed
##' @param ... Additional arguments passed to cor() or corKendall()
##' @return (d,d)-matrix of parameters *or* a d(d-1)/2 parameter vector
fitCor <- function(cop, x, method = c("kendall", "spearman"),
                   posDef = is(cop, "ellipCopula"), matrix = TRUE, ...) {
    method <- match.arg(method)
    theta <-
        switch(method,
               "kendall" = {
            tau <- corKendall(x, ...)
            iTau(cop, P2p(tau))
        },
        "spearman" = {
            rho <- cor(x, method=method, ...)
            iRho(cop, P2p(rho))
        },
        stop("Not yet implemented:", method))
    if (posDef) { # make pos. definite
        m <- nearPD(p2P(theta, d=ncol(x)), corr=TRUE)$mat # "dpoMatrix"
        if(!matrix) P2p(m) else m
    } else
        if(matrix) p2P(theta, d=ncol(x)) else theta
}

##' @title Computing the Log-Likelihood of the Given Copula
##' @param param Parameter (vector)
##' @param u The data in [0,1]^d
##' @param copula Copula object
##' @return Log-likelihood of the given copula at param given the data x
loglikCopula <- function(param, u, copula) {
    stopifnot((l <- length(param)) == length(copula@parameters)) # sanity check
    copula@parameters <- param
    lower <- copula@param.lowbnd
    if(length(lower) != l) return(-Inf)
    upper <- copula@param.upbnd
    if(length(upper) != l) return(-Inf)
    cop.param <- copula@parameters
    admissible <- !(any(is.na(cop.param) | cop.param > upper | cop.param < lower))
    if (admissible) {
        sum(dCopula(u, copula=copula, log=TRUE, checkPar=FALSE))
    } else -Inf
}

## TODO: Currently "unfinished";  instead  L <- .....(X)   where X <- getXmat() below
##' @title Variance-Covariance Matrix (vcov) for Inversion of Kendall's tau Estimator
##' @param copula The *fitted* copula object
##' @return L
getL <- function(copula) {
    ## for ellipCopula only!
    p <- copula@dimension
    pp <- p * (p - 1) / 2

    dgidx <- outer(1:p, 1:p, "-")
    dgidx <- P2p(dgidx)

    if (!is(copula, "ellipCopula") || copula@dispstr == "ex") {
        cbind(rep.int(1/pp, pp), deparse.level=0L)
    } else switch(copula@dispstr,
                  "un" = diag(pp),
                  "toep" = {
                     mat <- model.matrix(~ factor(dgidx) - 1)
                     mat / matrix(colSums(mat), nrow = pp, ncol=p - 1, byrow=TRUE)
                  },
                  "ar1" = {
               ## FIXME More efficient: estimate log(rho) first and then exponentiate back,
               ##  see e.g. fitCor(*, "spearman") above
               ## mat <- model.matrix(~ factor(dgidx) - 1)
               ## mat * matrix(1:(p - 1), nrow=pp, ncol=p - 1, byrow=TRUE)
               X <- getXmat(copula)
               ## L:
               t(solve(crossprod(X), t(X)))
           },
           stop("Not implemented yet for the dispersion structure ", copula@dispstr))
}

##' @title Design Matrix for Method-of-Moments Estimators via lm()
##' @param copula Copula object
##' @return Design matrix for method-of-moments estimators via lm()
getXmat <- function(copula) {
    p <- copula@dimension
    pp <- p * (p - 1) / 2
    if (!is(copula, "ellipCopula")) # one-parameter non-elliptical copula
	matrix(1, nrow=pp, ncol=1)
    else {
	switch(copula@dispstr,
	       "ex" = matrix(1, nrow=pp, ncol=1),
	       "un" = diag(pp),
	       "toep" =, "ar1" = {
            dgidx <- outer(1:p, 1:p, "-")
            dgidx <- P2p(dgidx)
            if(copula@dispstr == "toep")
                model.matrix(~ factor(dgidx) - 1)
            else { ## __"ar1"__
                ## estimate log(rho) first and then exponentiate back
                ## mat <- model.matrix(~ factor(dgidx) - 1)
                ## mat %*% diag(1:(p - 1))
                cbind(dgidx, deparse.level=0L)
            }
        },
        stop("Not implemented yet for the dispersion structure ", copula@dispstr))
    }
}

##' Auxiliary function for varKendall() and varSpearman()
##' @param cop An ellipCopula with \code{dispstr = "ar1"}
##' @param v Influence for tau or rho
##' @param L L
##' @param der A character string, currently either "tau" or "rho"
varInfluAr1 <- function(cop, v, L, der) {
    p <- cop@dimension
    pp <- p * (p - 1) / 2
    ## n <- nrow(v)

    ## estimate log(r) first, then log(theta), and then exponentiate back
    ## r is the lower.tri of sigma
    sigma <- getSigma(cop) # assuming cop is the fitted copula
    ## influ for log(r)
    r <- P2p(sigma)
    der <- if (der == "tau") dTauFun(cop)(r) else dRhoFun(cop)(r)
    D <- diag(x = 1 / r / der, pp)
    v <- v %*% D
    ## influ for log(theta)
    v <- v %*% L

    ## influ for theta
    theta <- cop@parameters[1]
    v %*% theta
}

## TODO this should *not* be used anymore; use Jscore() from gofCopula.R
## Auxiliary Function for varPL()
influ.terms <- function(u, influ, q)
{
    p <- ncol(u)
    n <- nrow(u)
    ## Integral w.r.t. empirical copula (-> sum):
    S <- matrix(0,n,q)
    for (i in 1:p) {
        o.i <- order(u[,i], decreasing=TRUE)
        obi <- rank(u[,i]) # FIXME? add.influ() in ./gofTests.R uses ecdf(.) * M here
        S <- S + rbind(rep.int(0,q),
                       apply(influ[[i]][o.i,,drop=FALSE],2,cumsum))[n + 1 - obi,,drop=FALSE]
    }
    S / n
}


### Method-of-moments based estimation #########################################

##' @title Variance of the Estimator based on Spearman's rho
##' @param cop The *fitted* copula object
##' @param u The data in [0,1]^d
##' @return Variance of the inversion of Spearman's rho estimator
##' @note See Kojadinovic & Yan (2010) "Comparison of three semiparametric ...",
##'       IME 47, 52--63
varSpearman <- function(cop,u)  {
    ## check if variance can be computed
    if (!hasMethod("dRho", class(cop))) {
        warning("The variance estimate cannot be computed for this copula.")
        q <- length(cop@parameters)
        return(matrix(NA, q, q))
    }

    stopifnot((p <- cop@dimension) >= 2)
    n <- nrow(u)
    v <- matrix(0, n, p*(p-1)/2)

    ord <- apply(u, 2, order, decreasing=TRUE)
    ordb <- apply(u, 2, rank)# ties : "average"
    storage.mode(ordb) <- "integer" # as used below

    l <- 0L
    for (j in 1:(p-1)) {
        for (i in (j+ 1L):p)
            v[,(l <- l + 1L)] <- u[,i] * u[,j] +
                c(0, cumsum(u[ord[,i], j]))[n + 1L - ordb[,i]] / n +
                c(0, cumsum(u[ord[,j], i]))[n + 1L - ordb[,j]] / n
    }
    ## X <- getXmat(cop)
    ## L <- t(solve(crossprod(X), t(X)))
    L <- getL(cop)
    v <- if (is(cop, "ellipCopula") && cop@dispstr == "ar1") {
        varInfluAr1(cop, v, L, "rho")
    } else {
        ## Caution: diag(0.5) is not a 1x1 matrix of 0.5!!! check it out.
        D <- if (length(cop@parameters) == 1) 1 / dRho(cop) else diag(1 / dRho(cop))
        v %*% L %*% D
    }
    144 * var(v)
}

##' @title Inversion of Spearman's rho estimator
##' @param copula The copula to be fitted
##' @param x The data in [0,1]^d or IR^d
##' @param estimate.variance A logical indicating whether the variance
##'        of the estimator shall be computed
##' @param warn.df A logical indicating whether a warning is given
##'        if the copula is coerced to df.fixed=TRUE
##' @param posDef A logical indicating whether a proper correlation matrix
##'        is computed
##' @param ... Additional arguments passed to fitCor()
##' @return The fitted copula object
fitCopula.irho <- function(copula, x, estimate.variance, warn.df=TRUE,
                           posDef=is(copula, "ellipCopula"), ...) {
    ccl <- getClass(class(copula))
    isEll <- extends(ccl, "ellipCopula")
    if(has.par.df(copula, ccl, isEll)) { ## must treat it as "df.fixed=TRUE"
        if(warn.df)
            warning("\"irho\" fitting ==> copula coerced to 'df.fixed=TRUE'")
        copula <- as.df.fixed(copula, classDef = ccl)
    }
    q <- length(copula@parameters)
    irho <- fitCor(copula, x, method = "spearman", posDef=posDef, matrix=FALSE, ...)
    X <- getXmat(copula)
    estimate <-
        as.vector(# stripping attributes
            if (isEll && copula@dispstr == "ar1") ## special treatment
                exp(.lm.fit(X, y=log(irho))$coefficients)
            else .lm.fit(X, y=irho)$coefficients)
    copula@parameters <- estimate
    var.est <- if (is.na(estimate.variance) || estimate.variance)
                   varSpearman(copula, x)/nrow(x)
               else matrix(NA, q, q)
    new("fitCopula",
        estimate = estimate,
        var.est = var.est,
        method = "inversion of Spearman's rho",
        loglik = NA_real_,
        fitting.stats = list(convergence = NA_integer_),
        nsample = nrow(x),
        copula = copula)
}

##' @title Variance of the Estimator based on Kendall's rho
##' @param cop the *fitted* copula object
##' @param u The data in [0,1]^d
##' @return variance of the inversion of Kendall's rho estimator
##' @note See Kojadinovic & Yan (2010) "Comparison of three semiparametric ...",
##'       IME 47, 52--63
varKendall <- function(cop, u) {
    ## check if variance can be computed
    if (!hasMethod("dTau", class(cop))) {
        warning("The variance estimate cannot be computed for a copula of class ", class(cop))
        q <- length(cop@parameters)
        return(matrix(NA, q, q))
    }
    p <- cop@dimension
    n <- nrow(u)
    ec <- numeric(n)
    v <- matrix(0,n,p*(p-1)/2)

    ## get influence functions for tau
    l <- 1
    for (j in 1:(p-1)) {
        for (i in (j+1):p) {
            for (k in 1:n) ## can this be vectorized?
                ec[k] <- sum(u[,i] <= u[k,i] & u[,j] <= u[k,j])/n
            v[,l] <- 2 * ec - u[,i] - u[,j]
            l <- l + 1
        }
    }
    ## X <- getXmat(cop)
    ## L <- t(solve(crossprod(X), t(X)))
    L <- getL(cop)
    v <- if (is(cop, "ellipCopula") && cop@dispstr == "ar1") { # special treatment
        varInfluAr1(cop, v, L, "tau")
    }
  else {
      ## Caution: diag(0.5) is not a 1x1 matrix of 0.5!!! check it out.
      D <- if (length(cop@parameters) == 1) 1 / dTau(cop) else diag(1 / dTau(cop))
      v %*% L %*% D
  }
    16 * var(v)
}

##' @title Inversion of Kendall's rho estimator
##' @param copula The copula to be fitted
##' @param x The data in [0,1]^d or IR^d
##' @param estimate.variance A logical indicating whether the variance
##'        of the estimator shall be computed
##' @param warn.df A logical indicating whether a warning is given
##'        if the copula is coerced to df.fixed=TRUE
##' @param posDef A logical indicating whether a proper correlation matrix
##'        is computed
##' @param ... Additional arguments passed to fitCor()
##' @return The fitted copula object
fitCopula.itau <- function(copula, x, estimate.variance, warn.df=TRUE,
                           posDef=is(copula, "ellipCopula"), ...) {
    ccl <- getClass(class(copula))
    isEll <- extends(ccl, "ellipCopula")
    if(has.par.df(copula, ccl, isEll)) { ## must treat it as "df.fixed=TRUE"
        ## currently: tCopula  or tevCopula
        if(warn.df)
            warning("\"itau\" fitting ==> copula coerced to 'df.fixed=TRUE'")
        copula <- as.df.fixed(copula, classDef = ccl)
    }
    q <- length(copula@parameters)
    itau <- fitCor(copula, x, method = "kendall", posDef=posDef, matrix=FALSE, ...)
    X <- getXmat(copula)
    estimate <-
        as.vector(# stripping attributes
            if(isEll && copula@dispstr == "ar1") ## special treatment
                exp(lm.fit(X, y=log(itau))$coefficients)
            else lm.fit(X, y=itau)$coefficients)
    copula@parameters <- estimate
    var.est <- if (is.na(estimate.variance) || estimate.variance)
                   varKendall(copula, x) / nrow(x) else matrix(NA, q, q)
    new("fitCopula",
        estimate = estimate,
        var.est = var.est,
        method = "inversion of Kendall's tau",
        loglik = NA_real_,
        fitting.stats = list(convergence = NA_integer_),
        nsample = nrow(x),
        copula = copula)
}


### "Inversion of Kendall's tau + MLE for nu" estimator for t copula ###########

##' @title Estimator of Demarta, McNeil (2005) for t Copulas
##' @param copula The copula to be fitted
##' @param u The data in [0,1]^d
##' @param posDef A logical indicating whether a proper correlation matrix
##'        is computed
##' @param lower The lower bound for optimize() (default 0 means Gaussian as we go in 1/nu)
##' @param upper The upper bound for optimize() (default 32 means down to nu = 1/32)
##' @param estimate.variance A logical indicating whether the estimator's
##'        variance shall be computed (TODO: not fully implemented yet)
##' @param hideWarnings A logical indicating whether warnings from optimize()
##'        shall be suppressed
##' @param tol Tolerance of optimize() for estimating nu
##' @param ... Additional arguments passed to fitCopula.itau()
##' @return The fitted copula object
##' @author Marius Hofert and Martin Maechler
fitCopula.itau.ml <- function(copula, u, posDef=TRUE, lower=NULL, upper=NULL,
                              estimate.variance, hideWarnings,
                              tol=.Machine$double.eps^0.25, ...)
{
    if(any(u < 0) || any(u > 1))
        stop("'u' must be in [0,1] -- probably rather use pobs(.)")
    if(!is(copula, "tCopula")) stop("method \"itau.ml\" is only applicable for \"tCopula\"")
    if(copula@df.fixed) stop("Use method=\"itau\" for 'tCopula' with 'df.fixed=TRUE'")
    stopifnot(is.numeric(d <- ncol(u)), d >= 2)
    if (copula@dimension != d)
        stop("The dimension of the data and copula do not match")
    if(is.null(lower)) lower <- 0  # <=>  df=Inf <=> Gaussian
    if(is.null(upper)) upper <- 32 # down to df = 1/32

    ## Estimation of the t-copula with the approach of Demarta, McNeil (2005)
    ## 1) Estimate the correlation matrix P
    fm <- fitCopula.itau(copula, u, estimate.variance=FALSE,
                         warn.df=FALSE, # (to change it there)
                         posDef=posDef, ...)
    p. <- fm@estimate
    ## 2) Estimate the d.o.f. parameter nu  via Maximum Likelihood :
    if(FALSE) # for debugging
        nLL <- function(Inu) {
            r <- -loglikCopula(c(p., df=1/Inu), u=u, copula=copula)
            cat(sprintf("1/nu=%12g, df=%8.3g, logL=%g\n", Inu, 1/Inu, -r))
            r
        }
    nLL <- function(Inu) -loglikCopula(c(p., df=1/Inu), u=u, copula=copula) # negative log-likelihood in 1/nu given P
    (if(hideWarnings) suppressWarnings else identity)(
        fit <- optimize(nLL, interval=c(lower, upper), tol=tol) # optimize
    )
    q <- length(copula@parameters) # the last is nu
    copula@parameters[seq_len(q-1L)] <- p.
    copula@parameters[[q]] <- copula@df <- df <- 1/fit$minimum
    loglik <- fit$objective
    has.conv <- TRUE # FIXME? use tryCatch() above to catch non-convergence
    if (is.na(estimate.variance))
        estimate.variance <- FALSE ## not yet: has.conv
    ## if(!has.conv)
    ##     warning("possible convergence problem: . . . . . . .)

    varNA <- matrix(NA_real_, q, q)
    var.est <- if(estimate.variance) {
        stop("'estimate.variance' not yet implemented for  \"itau.ml\" method")
        ## TODO: use one/zero step of fitCopula.mpl(*...., maxit=0) to get full vcov()
    } else varNA

    new("fitCopula",
        estimate = c(p., df),
        var.est = var.est,
        method = "itau for dispersion matrix P and maximum likelihood for df",
        loglik = loglik,
        fitting.stats = list(convergence = 0),
        ## optimize() does not give any info! -- if we had final optim(*, hessian=TRUE) step?
        ## c(list(method=method),
        ## fit[c("convergence", "counts", "message")], control),
        nsample = nrow(u),
        copula = copula)
}


### Maximum likelihood estimation ##############################################

##' @title Computing Initial Values for fitCopula.ml()
##' @param copula The copula to be fitted
##' @param u The data in [0,1]^d
##' @param default The default initial values
##' @param ... Additional arguments passed to fitCopula.itau()
##' @return Initial value for fitCopula.ml()
fitCopStart <- function(copula, u, default=copula@parameters, ...)
{
    clc <- class(copula)
    if(hasMethod("iTau", clc)) {
	ccl <- getClass(clc)
	.par.df <- has.par.df(copula, ccl)
	start <- fitCopula.itau(if(.par.df) as.df.fixed(copula, ccl) else copula,
				x=u, estimate.variance=FALSE, warn.df=FALSE, ...)@estimate # fitCopula.itau()
	if(.par.df) # add starting value for 'df'
	    start <- c(start, copula@df)
	if(is.finite(loglikCopula(start, u=u, copula=copula))) start else default
    } else default
}

##' @title Maximum Likelihood Estimator for Copulas
##' @param copula The copula to be fitted
##' @param u The data in [0,1]^d
##' @param start The initial value for optim()
##' @param lower The vector of lower bounds for optim()
##' @param upper The vector of upper bounds for optim()
##' @param optim.method The optimization method for optim()
##' @param optim.control optim()'s control parameter
##' @param estimate.variance A logical indicating whether the estimator's
##'        variance shall be computed
##' @param hideWarnings A logical indicating whether warnings from optim()
##'        shall be suppressed
##' @param bound.eps A small quantity denoting an eps for staying away from
##'        the theoretical parameter bounds
##' @param ... Additional arguments (currently with no effect)
##' @return The fitted copula object
fitCopula.ml <- function(copula, u, start, lower, upper,
                         optim.method, optim.control, estimate.variance, hideWarnings,
                         bound.eps=.Machine$double.eps^0.5, ...)
{
    chk.s(...)
    if(any(u < 0) || any(u > 1))
        stop("'u' must be in [0,1] -- probably rather use pobs(.)")
    stopifnot(is.numeric(d <- ncol(u)), d >= 2)
    if (copula@dimension != d)
        stop("The dimension of the data and copula do not match")
    if(is.null(start))
        start <- fitCopStart(copula, u=u)
    if(any(is.na(start))) stop("'start' contains NA values")
    q <- length(copula@parameters)
    if (q != length(start))
        stop(gettextf("The lengths of 'start' (= %d) and copula@parameters (=%d) differ",
                      length(start), q), domain=NA)

    control <- c(optim.control, fnscale = -1)
    control <- control[!vapply(control, is.null, NA)]

    if (!is.null(optim.control[[1]])) control <- c(control, optim.control)
    meth.has.bounds <- optim.method %in% c("Brent","L-BFGS-B")
    if (is.null(lower))
        lower <- if(meth.has.bounds) copula@param.lowbnd + bound.eps else -Inf
    if (is.null(upper))
        upper <- if(meth.has.bounds) copula@param.upbnd  - bound.eps else Inf

    ## messageOut may be used for debugging
    ## if (hideWarnings) {
    ## ## ?sink mentions that sink(*, type="message") is too dangerous: (!)
    ## As we found too, it swallows error messages and makes fitCop*() return nothing!

    ## messageOut <- textConnection("fitMessages", open="w", local=TRUE)
    ## sink(messageOut)
    ## sink(messageOut, type="message")
    ## oop <- options(warn = -1) ## ignore warnings; can be undesirable!
    ## on.exit({ options(oop); sink(type="message"); sink(); close(messageOut)})
    ## }

    (if(hideWarnings) suppressWarnings else identity)(
        fit <- optim(start, loglikCopula,
                     lower=lower, upper=upper,
                     method=optim.method,
                     copula = copula, u = u,
                     control = control)
    )

    ## if (hideWarnings) {
    ##   options(oop); sink(type="message"); sink()
    ##   on.exit()
    ##   close(messageOut)
    ## }

    copula@parameters[1:q] <- fit$par
    loglik <- fit$val
    has.conv <- fit[["convergence"]] == 0
    if (is.na(estimate.variance))
        estimate.variance <- has.conv
    if(!has.conv)
        warning("possible convergence problem: optim() gave code=",
                fit$convergence)

    ## MM{FIXME}: This should be done only by 'vcov()' and summary() !
    varNA <- matrix(NA_real_, q, q)
    var.est <- if(estimate.variance) {
        fit.last <- optim(copula@parameters, loglikCopula, lower=lower, upper=upper,
                          method=optim.method, copula=copula, u=u,
                          control=c(control, maxit=0), hessian=TRUE)
        vcov <- tryCatch(solve(-fit.last$hessian), error = function(e) e)
        if(is(vcov, "error")) {
            warning("Hessian matrix not invertible: ", vcov$message)
            varNA
        } else vcov ## ok
    } else varNA

    new("fitCopula",
        estimate = fit$par,
        var.est = var.est,
        method = "maximum likelihood",
        loglik = loglik,
        fitting.stats = c(list(method=optim.method),
                          fit[c("convergence", "counts", "message")], control),
        nsample = nrow(u),
        copula = copula)
}


### Maximum pseudo-likelihood estimation #######################################

##' @title Variance-Covariance Matrix (vcov) for Pseudo Likelihood Estimate
##' @param cop The *fitted* copula object
##' @param u The data in [0,1]^d
##' @return vcov matrix
varPL <- function(cop, u)
{
    q <- length(cop@parameters)
    p <- cop@dimension
    ans <- matrix(NA_real_, q, q)
    ccl <- getClass(clc <- class(cop))
    isEll <- extends(ccl, "ellipCopula")
    ## check if variance can be computed
    msg <- gettext("The variance estimate cannot be computed for this copula.",
                   " Rather use 'estimate.variance = FALSE'")
    if(is(cop, "archmCopula")) {
	fam <- names(which(.ac.classNames == class(cop)[[1]]))
	msg <- c(msg, gettext(" Or rather  emle(u, oCop)  instead; where",
                              sprintf(" oCop <- onacopula(%s, C(NA, 1:%d))", fam, p)))
    }
    if (!isEll && (!hasMethod("dcopwrap", clc) ||
                    !hasMethod("derPdfWrtArgs", clc))) {
	warning(msg); return(ans)
    }

    ## n <- nrow(u)

    ## influence: second part
    ## integrals computed from the original pseudo-obs u by Monte Carlo
    dcop <- dcopwrap(cop,u) ## wrapper: either dCopula() or '1' (for ellip.)
    ## New    dcop <- dCopula(u,cop) ## in some cases
    ## influ0 <- score (cop,u, dcop)
    ## derArg <- score2(cop,u, dcop)
    influ0 <- derPdfWrtParams(cop,u) / dcop # c. / c  = of dim.  n x q  {== cop<foo>@score()}
    if(is.na(influ0[1])) {
	warning(msg); return(ans)
    }
    derArg <- derPdfWrtArgs(cop,u)   / dcop #         = of dim.  n x p  {missing in copFoo -- TODO}

    ## TODO: use  array instead of list (and change influ.terms() accordingly)
    influ <- vector("list",p)
    for (i in 1:p)
        influ[[i]] <- influ0 * derArg[,i]

    inve <- solve(var(influ0))
    if(has.par.df(cop, ccl, isEll)) {
        ## currently cannot get var/cov for 'df' part
	ans[-q,-q] <- inve %*% var(influ0 - influ.terms(u,influ, q-1)) %*% inve
	ans
    }
    else
	inve %*% var(influ0 - influ.terms(u,influ, q)) %*% inve
}

##' @title Maximum Pseudo-Likelihood Estimator for Copulas
##' @param copula The copula to be fitted
##' @param u The data in [0,1]^d
##' @param start The initial value for optim()
##' @param lower The vector of lower bounds for optim()
##' @param upper The vector of upper bounds for optim()
##' @param optim.method The optimization method for optim()
##' @param optim.control optim()'s control parameter
##' @param estimate.variance A logical indicating whether the estimator's
##'        variance shall be computed
##' @param hideWarnings A logical indicating whether warnings from optim()
##'        shall be suppressed
##' @param ... Additional arguments passed to fitCopula.ml()
##' @return The fitted copula object
fitCopula.mpl <- function(copula, u, start, lower, upper,
                          optim.method, optim.control, estimate.variance,
                          hideWarnings, ...)
{
    fit <- fitCopula.ml(copula, u, start=start, lower=lower, upper=upper,
                        optim.method=optim.method, optim.control=optim.control,
                        estimate.variance=FALSE, hideWarnings=hideWarnings, ...)
    if (is.na(estimate.variance))
        estimate.variance <- fit@fitting.stats[["convergence"]] == 0
    var.est <- if(estimate.variance) {
                   varPL(fit@copula, u) / nrow(u)
               } else {
                   q <- length(copula@parameters)
                   matrix(NA, q, q)
               }
    new("fitCopula",
        estimate = fit@estimate,
        var.est = var.est,
        method = "maximum pseudo-likelihood",
        loglik = fit@loglik,
        fitting.stats = fit@fitting.stats,
        nsample = fit@nsample,
        copula = fit@copula)
}


### Wrapper ####################################################################

##' @title Main Fitting Wrapper
##' @param copula The copula to be fitted
##' @param data The data in [0,1]^d for "mpl", "ml", "itau.ml"; for "itau", "irho",
##'        it can be in [0,1]^d or IR^d
##' @param method The estimation method
##' @param posDef A logical indicating whether pairwise estimated correlation
##'        matrices are turned into proper correlation matrices (via nearPD())
##' @param start The initial value for optim()
##' @param lower The vector of lower bounds for optim()
##' @param upper The vector of upper bounds for optim()
##' @param optim.method The optimization method for optim()
##' @param optim.control optim()'s control parameter
##' @param estimate.variance A logical indicating whether the estimator's
##'        variance shall be computed
##' @param hideWarnings A logical indicating whether warnings from optim()
##'        shall be suppressed
##' @param ... Additional arguments passed to auxiliary functions
##' @return The fitted copula object
fitCopula <- function(copula, data,
                      method=c("mpl", "ml", "itau", "irho", "itau.ml"),
                      posDef=is(copula, "ellipCopula"),
                      start=NULL, lower=NULL, upper=NULL,
                      optim.method="BFGS", optim.control=list(maxit=1000),
                      estimate.variance=NA, hideWarnings=TRUE, ...)
{
    if(!is.matrix(data)) {
        warning("coercing 'data' to a matrix.")
        data <- as.matrix(data); stopifnot(is.matrix(data))
    }
    method <- match.arg(method)
    switch(method,
           "ml" =
               fitCopula.ml(copula, data, start=start, lower=lower, upper=upper,
                            optim.method=optim.method, optim.control=optim.control,
                            estimate.variance=estimate.variance,
                            hideWarnings=hideWarnings, ...),
           "mpl" = # Note that fitCopula.mpl() calls *.ml():
               fitCopula.mpl(copula, data, start=start, lower=lower, upper=upper,
                             optim.method=optim.method, optim.control=optim.control,
                             estimate.variance=estimate.variance,
                             hideWarnings=hideWarnings, ...),
           "itau.ml" = fitCopula.itau.ml(
               copula, data, posDef=posDef, lower=lower, upper=upper,
               estimate.variance=estimate.variance,
               hideWarnings=hideWarnings, ...), # <- may include 'tol' !
           "itau" = fitCopula.itau(copula, data,
                                   estimate.variance=estimate.variance, ...),
           "irho" = fitCopula.irho(copula, data,
                                   estimate.variance=estimate.variance, ...))
}
