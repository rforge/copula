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

## ==== Blomqvist's beta =======================================================

##' Compute the sample version of Blomqvist's beta for an Archimedean copula,
##' see, e.g., Schmid and Schmidt (2007) "Nonparametric inference on multivariate 
##' versions of Blomqvist’s beta and related measures of tail dependence"
##' @param u matrix of realizations following the copula
##' @param cop acopula to be estimated
##' @param scaling if FALSE then the scaling factors 2^(d-1)/(2^(d-1)-1) and
##'                2^(1-d) are omitted
##' @return sample version of multivariate Blomqvist beta
##' @author Marius Hofert
beta.hat <- function(u, cop, scaling = TRUE){
    stopifnot(all(0 <= u, u < 1))
    less.u <- u <= 0.5
    prod1 <- apply(less.u,1,all)
    prod2 <- apply(!less.u,1,all)
    sum.prod <- prod1 + prod2
    b <- mean(sum.prod)
    d <- ncol(u)
    if(scaling) 2^(d-1)/(2^(d-1)-1)*(b-2^(1-d)) else b
}

##' Compute the population version of Blomqvist's beta for an Archimedean copula
##' @param cop acopula to be estimated
##' @param theta copula parameter
##' @param d dimension
##' @param scaling if FALSE then the scaling factors 2^(d-1)/(2^(d-1)-1) and
##'                2^(1-d) are omitted
##' @return population version of multivariate Blomqvist beta 
##' @author Marius Hofert
beta. <- function(cop, theta, d, scaling = TRUE){
    j <- 1:d
    signs <- (-1)^j
    ldiags <- cop@diag(0.5,theta,j,log=TRUE)
    b <- cop@diag(0.5,theta,d) + 1 + sum(signs*exp(lchoose(d,j)+ldiags))
    if(scaling) 2^(d-1)/(2^(d-1)-1)*(b-2^(1-d)) else b
}

##' Method-of-moment-like estimation of Archimedean copulas based on a  
##' multivariate version of Blomqvist's beta
##' @param u matrix of realizations following the copula
##' @param cop acopula to be estimated
##' @param ... additional parameters for optimx()
##' @return Blomqvist beta estimator
##' @author Marius Hofert
ebeta <- function(u,cop,...){
    
    ## Note: We do not need the constants 2^(d-1)/(2^(d-1)-1) and 2^(1-d) here,
    ##	     since we equate the population and sample versions of Blomqvist's 
    ##       beta anyway.
    b.hat <- beta.hat(u, cop, scaling = FALSE)
    
    ## objective function
    d <- ncol(u)
    distance <- function(theta) abs(beta.(cop,theta,d,scaling = FALSE)-b.hat)

    ## optimize
    ## for the initial value, use DMLE for Gumbel and convert the parameters via tau
    start <- cop@tauInv(copGumbel@tau(edmle(u,copGumbel,...)))
    optimx::optimx(start,distance,lower = min(cop@paraInterval),upper = max(
                                                        cop@paraInterval),...)

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
etau <- function(u,cop,method = c("tau.mean","theta.mean"),...){
    
    tau.hat.mat <- cor(u,method="kendall",...) # matrix of pairwise tau()
    tau.hat <- tau.hat.mat[upper.tri(tau.hat.mat)] # all tau hat's
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

##' Maximum likelihood estimation based on the diagonal of the Archimedean copula
##' @param u matrix of realizations following the copula
##' @param cop acopula to be estimated
##' @param ... additional parameters for optimx()
##' @return diagonal maximum likelihood estimator
##' @author Marius Hofert
edmle <- function(u,cop,...){
    
    x <- apply(u,1,max) # data from the diagonal
    d <- ncol(u) # dimension
    l.d <- log(d) 
    ## explicit estimator for Gumbel (used as initial value for the other families as well)
    dmle.hat.G <- l.d/(log(length(x))-log(sum(-log(x)))) 
    if(cop@name == "Gumbel"){
	dmle.hat.G
    }else{
	tau.hat.G <- copGumbel@tau(dmle.hat.G)
	if(cop@name == "AMH" && tau.hat.G >= 0.33){ # AMH cannot attain a tau >= 1/3
            start <- 0.99 # most likely, this does not converge then (especially if tau.hat.G >> 0.33)
	}else{
            start <- cop@tauInv(tau.hat.G) # convert tau.hat.G to the parameter of cop
        }
        ## optimize
	mDiagLogL <- function(theta) sum(cop@dDiag(x,theta,d,log = TRUE)) # -log-Likelihood of the diagonal
	optimx::optimx(start,mDiagLogL,lower=min(cop@paraInterval),upper=max(
                                                           cop@paraInterval),...) 
    }

} 

## ==== (Simulated) maximum likelihood estimation ==============================

##' (Simulated) maximum likelihood estimation for Archimedean copulas
##' @param u matrix of realizations following the copula
##' @param cop acopula to be estimated
##' @param method estimation method
##' @param N approximation parameter for MLE; sample size for SMLE
##' @param initial method used for finding an initial value (either etau with 
##' 	   option "tau.mean" or "theta.mean" (denoted by "tau.mean" and "tau.mean", 
##'	   respectively), or edmle (denoted by "diag"))
##' @param ... additional parameters for optimx()
##' @return (simulated) maximum likelihood estimator
##' @author Marius Hofert
emle <- function(u,cop,method = c("mle","smle"),N,initial = c("tau.mean",
                                                  "theta.mean","diag"),...){
    ## compute initial value based on either etau or edmle
    start <- switch(initial,
                    tau.mean = {etau(u,cop,method="tau.mean")},
                    theta.mean = {etau(u,cop,method="theta.mean")},
                    diag = {edmle(u,cop,...)},
                {stop("wrong argument initial")})
    
    ## optimize
    mLogL <- function(theta) -sum(dnacopula(u,cop,theta,MC <- if(method == "smle") 
                                            TRUE else FALSE, N, log = TRUE)) # -log-Likelihood 
    optimx::optimx(start,mLogL,lower=min(cop@paraInterval),upper=max(
                                                   cop@paraInterval),...)
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
##'        "mle.tau.mean"     MLE with initial value found by averaged pairwise Kendall's tau       	
##'        "mle.theta.mean"   MLE with initial value found by Kendall's tau estimators averaged
##'        "mle.diag"         MLE with initial value found via DMLE
##'        "smle.tau.mean"    SMLE with initial value found by averaged pairwise Kendall's tau       
##'        "smle.theta.mean"  SMLE with initial value found by Kendall's tau estimators averaged
##'        "smle.diag"        SMLE with initial value found via DMLE
##'        "tau.tau.mean"     averaged pairwise Kendall's tau estimator 
##'        "tau.theta.mean"   average of Kendall's tau estimators
##'        "dmle"             MLE based on the diagonal
##'        "beta"             multivariate Blomqvist's beta estimator
##' @param ... additional parameters for optimx()
##' @return estimator according to the chosen method
##' @author Marius Hofert
enacopula <- function(x,cop,method = c("mle.tau.mean","mle.theta.mean","mle.diag",
                            "smle","tau.tau.mean","tau.theta.mean","dmle","beta"),
                      N,do.pseudo = FALSE,...){

    ## setup cop
    if(class(cop) != "outer_nacopula"){
	stop("cop has to be of class outer_nacopula")
    }
    if(cop@childCops != list()){
	stop("currently, only Archimedean copulas are provided")
    }

    ## setup x
    if(do.pseudo) u <- pobs(x) else u <- x

    ## main part
    switch(method,
           mle.tau.mean = {emle(u,cop,"mle",N,"tau.mean",...)},
           mle.theta.mean = {emle(u,cop,"mle",N,"theta.mean",...)},
           mle.diag = {emle(u,cop,"mle",N,"diag",...)},
           smle.tau.mean = {emle(u,cop,"smle",N,"tau.mean",...)},
           smle.theta.mean = {emle(u,cop,"smle",N,"theta.mean",...)},
           smle.diag = {emle(u,cop,"smle",N,"diag",...)},
           tau.tau.mean = {etau(u,cop,"tau.mean",...)},
           tau.theta.mean = {etau(u,cop,"theta.mean",...)},
           dmle = {edmle(u,cop,...)},
           beta = {ebeta(u,cop,...)},
       {stop("wrong estimation method")})

}	
