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
    d <- ncol(u)
    u. <- switch(method,
                 "log" = { pgamma(rowSums(-log(u)),shape=d) },
                 "normal" = { pchisq(rowSums(pnorm(u)^2),d) },
                 stop("wrong choice of method"))
    ## check U[0,1] of u.
    ad.test(u.)
}

##' Transforms vectors of random variates following the given nested Archimedean 
##' copula (with specified parameters) to U[0,1]^d vectors of random variates
##' @param x data matrix
##' @param theta parameter theta
##' @param do.pseudo boolean indicating whether to compute the pseudo-observations
##' @param MC whether or not K is evaluated by Monte Carlo simulation
##' @param N approximation parameter for series truncation; sample size for MC
##' @return matrix of supposedly U[0,1]^d realizations
##' @author Marius Hofert
gnacopulatrafo <- function(x, cop, do.pseudo = FALSE, MC = FALSE, N){
    stopifnot((d <- ncol(x)) >= 2)
    if(cop@childCops != list()) stop("currently, only Archimedean copulas are provided")
    u <- if(do.pseudo) pobs(x) else x
    psiI <- cop@psiInv(u,cop@theta)
    denom <- psiI[,1]
    res <- matrix(, nrow = nrow(u), ncol = d)
    ## compute first d-1 components of the transformation
    for(j in seq_len(d-1)) {
        num <- denom
        denom <- num + psiI[,j+1]
        res[,j] <- (num/denom)^j
    }
    ## compute dth component
    res[,d] <- if(cop@name == "Clayton"){
	cop@K(denom, cop@theta, d)
    }else{
        cop@K(denom, cop@theta, d, MC, N)
    }
    ## return transformed data 
    res
}

##' Conducts a goodness-of-fit test for the given H0 copula cop based on the 
##' data x
##' @param x data matrix
##' @param cop nacopula with specified H0 parameters
##' @param do.pseudo boolean indicating whether to compute the pseudo-observations
##' @param ad.test boolean indicating whether g01 is applied to the transformed 
##'        data or if just the transformed data is returned
##' @param method for Anderson-Darling, see g01
##' @param MC whether or not K is evaluated by Monte Carlo simulation
##' @param N approximation parameter for series truncation; sample size for MC
##' @param bootstrap whether or not a bootstrap is applied
##' @param M number of bootstrap replications
##' @param estimation.method estimation method, see enacopula
##' @param verbose if TRUE, the progress of the bootstrap is displayed
##' @return either a matrix of supposedly U[0,1]^d realizations is returned
##'         (the default) or the ad.test result
##' FIXME: return htest object...
##' @author Marius Hofert
gnacopula <- function(x, cop, do.pseudo = TRUE, ad.test = TRUE, 
                      method = c("log","normal"), MC = FALSE, N, bootstrap = TRUE, 
                      M = 1000, estimation.method = c("mle.tau.mean",
                                "mle.theta.mean","mle.diag","smle","tau.tau.mean",
                                "tau.theta.mean","dmle","beta"),verbose = TRUE){
    
    stopifnot((d <- ncol(x)) >= 2)
    if(cop@childCops != list()) stop("currently, only Archimedean copulas are provided")		
    u <- if(do.pseudo) pobs(x) else x
    if(bootstrap && !ad.test) stop("bootstrap requires ad.test == TRUE")

    if(bootstrap){ 
	
	## bootstrap setup
	n <- nrow(u)
	## although not recommended, bootstrap can be used with do.pseudo == FALSE
	## in which case the data should come from the copula directly (so at least
	## it should be between 0 and 1)
	stopifnot(all(0 <= u, u <= 1))
	## estimate the parameter by the provided method 
	theta.hat <- enacopula(u,cop,estimation.method,N,do.pseudo = do.pseudo) 
	## define the copula with theta.hat
	cop.hat <- onacopulaL(cop@family,list(theta.hat,1:d))
	## transform the data with the copula with estimated parameter and compute
	## the AD test result
	AD <- g01(gnacopulatrafo(u, cop.hat, do.pseudo = do.pseudo, MC = MC, N), 
                  method = method)

        ## conduct the parametric bootstrap
	AD.vec <- numeric(N) # vector of AD test statistics
	for(k in 1:N){
            
            ## do work
            u. <- rnacopula(n,cop.hat) # sample the copula with estimated parameter
            u.boot <- if(do.pseudo) pobs(u.) else u. # compute pseudo-observations if do.pseudo == TRUE
            theta.hat.boot <- enacopula(u.boot,cop,estimation.method,N,
                                        do.pseudo = do.pseudo)  # estimate the copula parameter
            cop.hat.boot <- onacopulaL(cop@family,list(theta.hat.boot,1:d)) # copula with estimated parameter
            AD.vec[k] <- g01(gnacopulatrafo(u.boot, cop.hat.boot, do.pseudo = 
                                            do.pseudo, MC = MC, N), method = method)

            ## progress output
            if(verbose) cat("bootstrap progress:",format(round(k*100/N),nsmall = 3,
                                                         trim = TRUE,
                                                         scientific = FALSE)," %\n",
                            sep="")
            
	}
	## estimate p-value
	mean(AD.vec > AD)
	
    }else{ # bootstrap == FALSE
	gnacopulatrafo(u, cop, do.pseudo = FALSE, MC = MC, N = N)
    }
    
}

