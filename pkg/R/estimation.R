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

## ==== Pairwise tau estimation ================================================

##' Averaged pairwise Kendall's tau estimator
##' @param u matrix of realizations following the copula
##' @param cop acopula to be estimated
##' @return averaged pairwise Kendall's tau estimator
##' @author Marius Hofert
etau <- function(u,cop){
    cor.mat <- cor(u,method="kendall") # matrix of pairwise Kendall's taus
    tau.hat.mean <- mean(cor.mat[upper.tri(cor.mat)]) # mean of the d*(d-1)/2 pairs
    cop@tauInv(tau.hat.mean)
}   

## ==== Diagonal maximum likelihood estimation =================================

##' Diagonal maximum likelihood estimator
##' @param u matrix of realizations following the copula
##' @param cop acopula to be estimated
##' @param ... additional parameters to optimize()
##' @return diagonal maximum likelihood estimator
##' @author Marius Hofert
edmle <- function(u,cop,...){
    d <- ncol(u)
    l.d <- log(d)
    x <- apply(u,1,max) # data from the diagonal
    if(cop@name == "Gumbel"){
        l.d/(log(length(x))-log(sum(-log(x)))) # estimator directly known
    }else{	
	f <- function(theta){
            psiI <- cop@psiInv(x,theta)
            l.d+cop@psiDAbs(d*psiI,theta,log = TRUE)-cop@psiDAbs(psiI,theta,log = TRUE)
	}
	interval <- c(0,.Machine$double.xmax) ## FIXME: should depend on the copula family
	optimize(-f,interval,...)$minimum	
    }
}


## ==== Simulated maximum likelihood estimation ================================

## ==== Maximum likelihood estimation ==========================================


##' Computes different parameter estimators for an Archimedean copula
##' @param u matrix of realizations following the copula
##' @param cop acopula to be estimated
##' @param method estimation method
##' @return estimator
##' @author Marius Hofert
eacopula <- function(u,cop,method="MLE"){
    ## remark: this is not optimal yet due to the following reasons:
    ## log(density()) evaluates the setup steps "not depending on theta" 
    ## 	   for every call---this is inefficient

    ## check data

    interval <- switch(cop@name,
                       AMH = {c(0.01,0.99)}, # upper bound corresponds to tau = 0.33
                       Clayton = {c(0.01,4.666667)}, # upper bound corresponds to tau = 0.7
                       Frank = {c(0.01,11.41157)}, # upper bound corresponds to tau = 0.7
                       Gumbel = {c(1.01,3.333333)}, # upper bound corresponds to tau = 0.7
                       Joe = {c(1.01,5.463749)}, # upper bound corresponds to tau = 0.7
                   {stop("wrong family")})

    optimize(-cop@logDensity,interval,tol=0.001)$minimum
}	
