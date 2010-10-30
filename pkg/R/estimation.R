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

#### Estimation for Archimedean copulas

## ==== Blomqvist's beta =======================================================



## TODO: implement diagonals of the families (C.diag) with option log
betae <- function(u,cop,...){

    h.d <- 2^(d-1)/(2^(d-1)-1)  
    h.b <- 2^(1-d) 

    ## sample version of Blomqvist's beta, see, e.g., Schmid and Schmidt (2007) 
    ## "Nonparametric inference on multivariate versions of Blomqvist’s beta and 
    ## related measures of tail dependence"
    n <- nrow(u)
    d <- ncol(u)
    less.u <- u <= 0.5
    prod1 <- apply(less.u,1,all)
    prod2 <- apply(!less.u,1,all)
    sum.prod <- prod1 + prod2
    mean.sum <- mean(sum.prod)
    if(mean.sum < h.b) mean.sum <- h.b # adjust in order to not getting negative beta.hat's
    beta.hat <- h.d * (mean.sum - h.b)

## set up objective function
distance <- function() abs(beta.hat-blomqvist.beta(d,family,theta,h.d,h.b))

    ## optimize
    optimize(blomqvist.distance,interval(family),u=u,family=family,h.d=h.d,
             h.b=h.b,beta.hat=beta.hat,...)$minimum

}

## ==== Diagonal maximum likelihood estimation =================================

##' Diagonal maximum likelihood estimator
##' @param u matrix of realizations following the copula
##' @param cop acopula to be estimated
##' @param ... additional parameters for optimx()
##' @return diagonal maximum likelihood estimator
##' @author Marius Hofert
dmle <- function(u,cop,...){
    d <- ncol(u)
    l.d <- log(d)
    x <- apply(u,1,max) # data from the diagonal
    dmle.hat.G <- l.d/(log(length(x))-log(sum(-log(x)))) # estimator for Gumbel
    if(cop@name == "Gumbel"){
	dmle.hat.G
    }else{
	tau.hat.G <- copGumbel@tau(dmle.hat.G)
	if(cop@name == "AMH" && tau.hat.G >= 0.33){ # AMH cannot attain a tau >= 1/3
            start <- 0.99
	}else{
            start <- cop@tauInv(tau.hat.G)
        }
	l.d <- log(d)
	diag.logL <- function(theta){
            psiI <- cop@psiInv(x,theta)
            l.d+cop@psiDAbs(d*psiI,theta,log=TRUE)-cop@psiDAbs(psiI,theta,log=TRUE)
	}
	optimx(start,-diag.logL,lower=min(cop@paraInterval),
               upper=max(cop@paraInterval),...) 
    }
}

## ==== Pairwise tau estimation ================================================

##' Compute matrix of averaged pairwise cor() estimators 
##' @param u matrix of realizations following the copula
##' @param cop nacopula to be estimated
##' @param ... additional arguments to cor()
##' @return averaged pairwise cor() estimators
##' @author Marius Hofert
taue <- function(u,cop,...){
    tau.n.mat <- cor(u,method="kendall",...) # matrix of pairwise tau()
    tau.n.mean <- mean(tau.n.mat[upper.tri(tau.n.mat)]) # mean of the d*(d-1)/2 pairs
    cop@tauInv(tau.n.mean)
}  

## TODO taue. fuer mean von schaetzern
## TODO das gleiche noch fuer pairwise MLE!
## TODO blomqvist implementen

## TODO: if everything works for ACs, include tools for nested Archimedean copula 
## TODO: include an optional argument which finds the "nearest" matrix of pairwise 
## estimators that fulfills the nesting condition (this can then be used as a starting 
## point for estimation methods for nacopulas)


## ==== (Simulated) maximum likelihood estimation ==============================

## todo: auch in gof noch ueberall Variablen "klein schreiben"
## nimm kendall's tau als startwerte...
# mle <- function(u,cop,method="mle"){
#     stopifnot()
# }

## ==== Estimation wrapper =====================================================

##' Computes different parameter estimators for an Archimedean copula
##' @param x matrix of realizations following the copula
##' @param cop acopula to be estimated
##' @param method estimation method
##' @return estimator
##' @author Marius Hofert
eacopula <- function(x,cop,method="MLE",do.pseudo=FALSE){
	# stopifnot(check that cop is acopula, check x, apply pobs...)
    ## remark: this is not optimal yet due to the following reasons:
    ## log(density()) evaluates the setup steps "not depending on theta" 
    ## 	   for every call---this is inefficient

    
    # optimize(-cop@logDensity,interval,tol=0.001)$minimum
}	
