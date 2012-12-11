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


### Various goodness-of-fit tests ##############################################


### Test statistics ############################################################

##' Test statistics for various goodness-of-fit tests of (supposedly) U[0,1]^d
##' distributed vectors of random variates
##'
##' @title Test statistics for U[0,1]^d
##' @param u matrix of supposedly  U[0,1]^d observations
##' @param method various test statistics. Available are:
##'        "Sn"     : the test statistic S_n (Cramer-von Mises) in Genest, Remillard, Beaudoin (2009)
##'        "SnB"    : the test statistic S_n^{(B)} in Genest, Remillard, Beaudoin (2009)
##'        "SnC"    : the test statistic S_n^{(C)} in Genest, Remillard, Beaudoin (2009)
##'        "AnChisq": Anderson-Darling test statistic after map to a chi-square
##'                   distribution using the standard normal quantile function
##'        "AnGamma": Anderson-Darling test statistic after map to an Erlang/Gamma
##'                   distribution using the logarithm
##' @param ... additional arguments for the different methods
##' @return values of the chosen test statistic
##' @author Marius Hofert and Martin Maechler
gofTstat <- function(u, method=c("Sn", "SnB", "SnC", "AnChisq", "AnGamma"), ...)
{
    if(!is.matrix(u)) u <- rbind(u, deparse.level=0L)
    d <- ncol(u)
    n <- nrow(u)
    method <- match.arg(method)
    switch(method,
           "Sn" =
       { ## S_n
           if(!hasArg(copula)) stop("object 'copula' required for pCopula() call")
           .C(cramer_vonMises,
              as.integer(n),
              as.integer(d),
              as.double(u),
              as.double(pCopula(u, ...)),
              stat=double(1))$stat
       },
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


### Computing different goodness-of-fit tests ##################################

##' Goodness-of-fit tests based on the parametric bootstrap
##'
##' @title Goodness-of-fit tests based on the parametric bootstrap
##' @param copula object of type 'copula' representing the H_0 copula
##'        (if necessary, parameters will be used as starting values for fitCopula())
##' @param x (n, d)-matrix containing the data
##' @param N number of bootstrap replications
##' @param method goodness-of-fit test statistic to be used; see ?gofTstat
##' @param estim.method estimation method for the unknown parameter vector; see ?fitCopula
##' @param optim.method optim() method used for fitting
##' @param optim.control see ?optim
##' @param trafo.method transformation to U[0,1]^d
##' @param verbose logical indicating whether a progress bar is shown
##' @param ... additional arguments passed to the transformation method
##' @return an object of class 'htest'
##' @author Ivan Kojadinovic, Marius Hofert
##' Note: - different '...' should be passed to different *.method
gofPB <- function(copula, x, N, method=eval(formals(gofTstat)$method),
                  estim.method=eval(formals(fitCopula)$method),
                  optim.method="BFGS", optim.control,
                  trafo.method=c("rtrafo", "htrafo"),
                  verbose=TRUE, ...)
{
    ## checks
    stopifnot(is(copula, "copula"), N>=1)
    if(!is.matrix(x)) x <- rbind(x, deparse.level=0L)
    stopifnot((d <- ncol(x))>1, (n <- nrow(x))>0, dim(copula)==d)
    method <- match.arg(method)
    estim.method <- match.arg(estim.method)
    trafo.method <- match.arg(trafo.method)
    if(trafo.method=="htrafo"){
        if(!is(copula, "outer_nacopula")) stop("trafo.method='htrafo' only implemented for copula objects of type 'outer_nacopula'")
        if(length(copula@childCops)) stop("currently, only Archimedean copulas are provided")
    }

    ## 1) compute the pseudo-observations
    uhat <- pobs(x)

    ## 2) fit the copula
    cop. <- fitCopula(copula, uhat, method=estim.method, estimate.variance=FALSE,
                      optim.method=optim.method, optim.control=optim.control)@copula

    ## 3) compute the test statistic
    u <- if(method=="Sn"){
        switch(trafo.method,
               "rtrafo"={
                   rtrafo(uhat, cop., ...)
               },
               "htrafo"={
                   htrafo(uhat, cop., ...)
               },
               stop("wrong transformation method"))
    } else uhat
    T <- if(method=="Sn") gofTstat(u, method=method, copula=cop.)
         else gofTstat(u, method=method)

    ## 4) simulation of the test statistic under the null hypothesis
    if(verbose){
	pb <- txtProgressBar(max=N, style=if(isatty(stdout())) 3 else 1) # setup progress bar
	on.exit(close(pb)) # on exit, close progress bar
    }
    T0 <- vapply(1:N, function(k){
        ## 4.1) sample the fitted copula
        Uhat <- pobs( rCopula(n, cop.) )

        ## 4.2) fit the copula
        cop.. <- fitCopula(copula, Uhat, method=estim.method, estimate.variance=FALSE,
			   optim.method=optim.method, optim.control=optim.control)@copula

        ## 4.3) compute the test statistic
        u. <- if(method=="Sn"){
            switch(trafo.method,
                   "rtrafo"={
                       rtrafo(Uhat, cop.., ...)
                   },
                   "htrafo"={
                       htrafo(Uhat, cop.., ...)
                   },
                   stop("wrong transformation method"))
        } else Uhat
        T0. <- if(method=="Sn") gofTstat(u., method=method, copula=cop..)
               else gofTstat(u., method=method)

        if(verbose) setTxtProgressBar(pb, k) # update progress bar
        T0. # return
    }, NA_real_)

    ## return result object
    structure(class = "htest",
	      list(method = sprintf(
                   "Parametric bootstrap goodness-of-fit test with 'method'=\"%s\", 'estim.method'=\"%s\"",
		   method, estim.method),
                   parameter = c(parameter = cop.@parameters),
                   statistic = c(statistic = T),
                   p.value = (sum(T0 >= T) + 0.5) / (N + 1), # typical correction =>  p-values in (0, 1)
                   data.name = deparse(substitute(x))))
}

##' Influence functions for the multiplier bootstrap
##'
##' @title Influence functions for the multiplier bootstrap
##' @param cop object of type 'copula'
##' @param u (n, d)-matrix of (pseudo-)observations
##' @param method "mpl" or one of "itau", "irho"
##' @return TODO: + rename
##' @author Marius Hofert (based on ideas of Ivan Kojadinovic)
##' Note: - References:
##'         * For estim.method="mpl":
##'           I. Kojadinovic and J. Yan (2011), A goodness-of-fit test for multivariate
##'           multiparameter copulas based on multiplier central limit theorems, Statistics
##'           and Computing 21:1, pages 17-30
##'         * For estim.method="itau" or ="irho" (d=2):
##'           I. Kojadinovic, J. Yan and M. Holmes (2011), Fast large-sample goodness-of-fit
##'           tests for copulas, Statistica Sinica 21:2, pages 841-871
Jscore <- function(cop, u, method)
{
    ## checks
    stopifnot(is(copula, "copula"))
    if(!is.matrix(u)) u <- rbind(u, deparse.level=0L)
    stopifnot((d <- ncol(x))>1, (n <- nrow(x))>0, dim(copula)==d)

    ## deal with different methods
    switch(method,
           "mpl"={ ## see page 22 in Kojadinovic and Yan (2011)
               ## integrals computed from n realizations by Monte Carlo
               dcop <- dcopwrap(cop, u)
               influ0 <- derPdfWrtParams(cop, u) / dcop
               derArg <- derPdfWrtArgs  (cop, u) / dcop
               influ <- lapply(1:cop@dimension, function(i) influ0*derArg[,i]) # TODO: improve
               n <- nrow(u)
               d <- ncol(u)
               p <- length(cop@parameters)
               S <- matrix(0, n, p)
               for(j in 1:d){ # TODO: remove for() and apply()
                   ij <- order(u[,j], decreasing=TRUE)
                   ijb <- vapply(1:n, function(i) sum(t(u[,j]) <= u[i,j]), NA_real_)
                   S <- S + rbind(rep.int(0, p),
                                  apply(influ[[j]][ij,,drop=FALSE], 2, cumsum))[n+1-ijb,,drop=FALSE] / n -
                                      matrix(colMeans(influ[[j]]*u[,j]), n, p, byrow=TRUE)
               }
               solve(crossprod(influ0)/n, t(derPdfWrtParams(cop, u)/dcopwrap(cop, u)-S))
           },
           "ml"={
               stop("method='ml' not available")
           },
           "itau"={ ## see page 849 in Kojadinovic, Yan, and Holmes (2011)
               stopifnot(dim(cop)==2)
               (4/dTau(cop)) * ( 2*pCopula(u, cop) - rowSums(u) + (1-tau(cop))/2 )
           },
           "irho"={ ## see page 847 in Kojadinovic, Yan, and Holmes (2011)
               stopifnot(dim(cop)==2)
               i1 <- order(u[,1], decreasing=TRUE)
               i2 <- order(u[,2], decreasing=TRUE)
               n <- nrow(u)
               i1b <- vapply(1:n, function(k) sum(t(u[,1]) <= u[k,1]), NA_real_)
               i2b <- vapply(1:n, function(k) sum(t(u[,2]) <= u[k,2]), NA_real_)
               rmu <- rowMeans(u)
               term <- c(0, cumsum(u[,2][i1]))[n+1-i1b] / n - rmu +
                   c(0, cumsum(u[,1][i2]))[n+1-i2b] / n - rmu
               (1/dRho(cop)) * ( 12*(apply(u, 1, prod) + term) - 3 - rho(cop) )
           },
           stop("wrong method"))
}

##' Goodness-of-fit tests based on the multiplier bootstrap
##'
##' @title Goodness-of-fit tests based on the multiplier bootstrap
##' @param copula object of type 'copula' representing the H_0 copula
##' @param x (n, d)-matrix containing the data
##' @param N number of bootstrap replications
##' @param method goodness-of-fit test statistic to be used; see ?gofTstat
##' @param estim.method estimation method for the unknown parameter vector; see ?fitCopula
##' @param optim.method optim() method used for fitting
##' @param optim.control see ?optim
##' @param ... additional arguments
##' @return an object of class 'htest'
##' @author Ivan Kojadinovic, Marius Hofert
##' Note: - Since we use C here (necessary?), we can't use verbose=TRUE
##'       - Also, since we use C, we can't easily pass different methods
##'         to gofTstat() => TODO
##'       - Apart from that, gofMP() is quite similar to gofPB()... (one could
##'         think about putting both in gofCopula()
gofMB <- function(copula, x, N, method=c("Sn"),
                  estim.method=eval(formals(fitCopula)$method),
                  optim.method="BFGS", optim.control, ...)
{
    ## checks
    stopifnot(is(copula, "copula"), N>=1)
    if(!is.matrix(x)) x <- rbind(x, deparse.level=0L)
    stopifnot((d <- ncol(x))>1, (n <- nrow(x))>0, dim(copula)==d)
    method <- match.arg(method)
    estim.method <- match.arg(estim.method)
    if(method!="Sn") stop("Multiplier method for '", method,"' not yet implemented")
    if(estim.method == "ml") stop("estim.method='ml' not available")
    if(estim.method %in% c("irho", "itau") && d > 2) stop("only bivariate case possible for estim.method='irho' or ='itau'")

    ## 1) compute the pseudo-observations
    uhat <- pobs(x)

    ## 2) fit the copula
    cop. <- fitCopula(copula, uhat, method=estim.method, estimate.variance=FALSE,
                      optim.method=optim.method, optim.control=optim.control)@copula

    ## 3) compute the test statistic
    T <- .C(cramer_vonMises_grid,
            as.integer(d),
            as.double(uhat),
            as.integer(n),
            as.double(uhat),
            as.integer(n),
            as.double(pCopula(uhat, cop.)),
            stat = double(1))$stat

    ## 4) simulation of the test statistic under H_0
    T0 <- .C(multiplier,
             as.integer(d),
             as.double(uhat),
             as.integer(n),
             as.double(uhat),
             as.integer(n),
             as.double(dCdtheta(cop., uhat) %*% Jscore(cop., u=uhat, method=estim.method)),
             as.integer(N),
             s0 = double(N))$s0

    ## return result object
    structure(class = "htest",
	      list(method = sprintf(
		   "Multiplier bootstrap goodness-of-fit test with 'method'=\"%s\", 'estim.method'=\"%s\"",
		   method, estim.method),
                   parameter = c(parameter = cop.@parameters),
                   statistic = c(statistic = T),
                   p.value = (sum(T0 >= T) + 0.5) / (N + 1), # typical correction =>  p-values in (0, 1)
                   data.name = deparse(substitute(x))))
}


### Wrapper ####################################################################

##' Goodness-of-fit test wrapper function
##'
##' @title Goodness-of-fit test wrapper function
##' @param copula object of type 'copula' representing the H_0 copula
##' @param x (n, d)-matrix containing the data
##' @param N the number of bootstrap (parametric or multiplier) replications
##' @param method goodness-of-fit test statistic to be used
##' @param estim.method estimation method for the unknown parameter vector
##' @param simulation parametric bootstrap ('pb') or multiplier method ('mult')
##' @param verbose logical indicating whether a progress bar is shown
##' @param print.every deprecated
##' @param optim.method optim() method used for fitting
##' @param optim.control see ?optim
##' @param ... additional arguments passed to the internal auxiliary functions
##'        gofPB() and gofMB()
##' @return an object of class 'htest'
##' @author Ivan Kojadinovic, Marius Hofert
gofCopula <- function(copula, x, N=1000,
                      method=eval(formals(gofTstat)$method),
                      estim.method=eval(formals(fitCopula)$method),
                      simulation=c("pb", "mult"),
		      verbose=TRUE, print.every=NULL,
                      optim.method="BFGS", optim.control=list(maxit=20), ...)
{
    ## checks
    stopifnot(is(copula, "copula"), N>=1)
    if(!is.matrix(x)) x <- rbind(x, deparse.level=0L)
    stopifnot((d <- ncol(x))>1, (n <- nrow(x))>0, dim(copula)==d)
    method <- match.arg(method)
    estim.method <- match.arg(estim.method)
    optim.method <- match.arg(optim.method)
    stopifnot(optim.method %in% eval(formals(optim)$method))
    ## deprecation
    if (!is.null(print.every)) {
        warning("Argument 'print.every' is deprecated. Please use 'verbose' instead.")
        verbose <- print.every > 0
    }
    ## back-compatibility
    if(missing(estim.method) && !missing(method)) {
        eMeth <- eval(formals()$estim.method)
	if(!is.na(i <- pmatch(method, eMeth))) {
	    warning("old (pre 0.999-*) argument 'method' is now called 'estim.method'")
	    estim.method <- eMeth[i]
	    method <- "Sn"
	}
    }

    ## distinguish the methods
    method <- match.arg(method)
    estim.method <- match.arg(estim.method)
    simulation <- match.arg(simulation)
    switch(simulation,
           "pb" = { ## parametric bootstrap
               gofPB(copula, x, N=N, method=method, estim.method=estim.method,
                     verbose=verbose, optim.method=optim.method,
                     optim.control=optim.control, ...)
           },
           "mult" = { ## multiplier bootstrap
               gofMB(copula, x=x, N=N, method=method, estim.method=estim.method,
                     optim.method=optim.method, optim.control=optim.control, ...)
           },
           stop("Invalid simulation method ", match.arg(simulation)))
}
