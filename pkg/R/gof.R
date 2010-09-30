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

#### Goodness-of-fit testing for (nested) Archimedean copulas

##' Computes the pseudo-observations for the given data matrix
##' @param X matrix of random variates to be converted to pseudo-observations
##' @return pseudo-observations (matrix of the same dimensions as X)
##' @author Marius Hofert
pobs <- function(X){
    n <- nrow(X)
    apply(X,2,function(x) rank(x)/(n+1))
}

##' Transforms U[0,1]^d vectors of random variates to U[0,1] distributed random
##' variates
##' @param U matrix of random variates to be transformed
##' @param method either "log" (map to an Erlang distribution) or "normal" (map
##' 	to a chi-square distribution)
##' @return U[0,1] distributed variates (the Erlang or chi-square cdf is applied)
##' @author Marius Hofert
tgofU <- function(U,method="log"){
    d <- ncol(U)
    switch(method,
           "log" = pgamma(rowSums(-log(U)),shape=d),
           "normal" = pchisq(rowSums(pnorm(U)^2),d),
           stop("wrong choice of method for tgofU()"))
}

##' Transforms d-dimensional vectors of random variates following the given
##' Archimedean copula (with specified parameter) to U[0,1]^d vectors of random
##' variates
##' @param X matrix of random variates 
##' @param cop acopula with specified parameter theta
##' @param copulaData boolean indicating whether X comes from a copula directly
##' @return matrix of the same dimension as X
##' @author Marius Hofert
tacopula <- function(X,cop,copulaData=FALSE){
    stopifnot(length(d <- dim(X)) == 2, (d <- d[2]) >= 2)
    if(!copulaData) U <- pobs(X) else U <- X
    theta <- cop@theta
    psiI.U <- cop@psiInv(U,theta)
    U.prime <- U # matrix(,nrow=n,ncol=d)
    sum.psiI.U <- psiI.U[,1]
    ## compute first d-1 components of the transformation
    for(j in seq_len(d-1)) {
        sum.. <- sum.psiI.U
        sum.psiI.U <- sum.. + psiI.U[,j+1]
        U.prime[,j] <- (sum../sum.psiI.U)^j
    }
    ## compute dth component
    U.prime[,d] <- cop@K(sum.psiI.U, theta, d)
    ## return transformed data
    U.prime
}

##' Conducts a goodness-of-fit test for the given H0 copula cop based on the 
##' data X
##' @title Goodness-Of-Fit Test for a Given (H0) Copula and Data
##' @param X matrix of data (for copulaData=TRUE, X should be data from a copula)
##' @param cop acopula with specified H0 parameter theta
##' @param copulaData boolean indicating whether X comes from a copula directly
##' @param N number of bootstrap replications
##' @param method applied for the mapping to one-dimensional data; either "log"
##'	(map to an Erlang distribution) or "normal" (map to a chi-square distribution)
##' @return p-value -- MM: FIXME: Rather "htest" object
##' @author Marius Hofert
gofacopula <- function(X,cop,copulaData=FALSE,N=1000,method="log", verbose=TRUE)
{
    if(copulaData){ # if data comes from a copula directly (no bootstrap is used)
        U <- tacopula(X,cop,copulaData) # transform the data to U[0,1]^d data under H0
        Y <- tgofU(U,method) # transform the U[0,1]^d data under H0 to U[0,1] data
        ad.test(Y)$p.value # compute the Anderson-Darling test statistic
    }else{ # pseudo-observations are computed and bootstrap is used
	## setup for bootstrap
	n <- nrow(X)
	d <- ncol(X)
	U. <- pobs(X) # compute pseudo-observations
	theta.hat <- mleacopula(U.,cop) # estimate the copula parameter
	cop.hat <- onacopulaL(cop@family,list(theta.hat,1:d)) # copula with estimated parameter
	U.. <- tacopula(U.,cop.hat,copulaData) # transform the data with the estimated parameter
	Y <- tgofU(U..,method) # map to one-dimensional data
	AD.stat <- ad.test(Y)$statistic # compute the Anderson-Darling test statistic
	## conduct the parametric bootstrap
	if(verbose) {
            cat(sprintf("parametric bootstrap (N = %d):\n",N))
            N. <- N %/% 10
        }
	AD.vec <- numeric(N)
	for(k in 1:N){
            ## progress output
            if(verbose & k%%N. == 1)
                cat(sprintf("bootstrap progress: %3.0f%%", k*100/N))
            ## do work
            U <- rnacopula(n,cop.hat) # sample the copula with estimated parameter
            U.star <- pobs(U) # compute pseudo-observations
            theta.hat.star <- mleacopula(U.star,cop) # estimate the copula parameter
            cop.hat.star <- onacopulaL(cop@family,list(theta.hat.star,1:d)) # copula with estimated parameter
            U..star <- tacopula(U.star,cop.hat.star,copulaData) # transform the data with the estimated parameter
            Y.star <- tgofU(U..star,method) # map to one-dimensional data
            AD.vec[k] <- ad.test(Y.star)$statistic # compute the Anderson-Darling test statistic
            if(verbose) ## progress output
                cat(",", if(k%%N. == 0) "\n")
	}
	## estimate p-value
	mean(AD.vec > AD.stat)
    }
}

# ##' Recursively transforms d-dimensional vectors of random variates following the given
# ##' nested Archimedean copula (with specified parameters) to U[0,1]^d vectors
# ##' of random variates
# ##' @param X matrix of random variates
# ##' @param cop nacopula with specified parameters
# ##' @param copulaData boolean indicating whether X comes from a copula directly
# ##' @return matrix of the same dimension as X
# ##' @author Marius Hofert
# tnacopula <- function(X,cop,copulaData=FALSE){
#     TODO
# First trial:
# if(!copulaData) U <- pobs(X) else U <- X
#     for all children do {
#         call procedure again
#     }
# consider parent

# }
