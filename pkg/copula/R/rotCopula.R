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

##################################################################################
### Rotated copulas created from an existing copula and a mask of logicals
##################################################################################

setClass("rotCopula", contains = "copula",
         representation = representation(
             copula = "copula",
             flip = "logical"
         ),
	 validity =
	     function(object) {
		 d <- object@copula@dimension
		 if(length(object@flip) != d)
		     "The dimension of the copula does not match the length of 'flip'"
		 else if(anyNA(object@flip))
		     "'flip' contains NA/NaN"
		 else TRUE
	     })

##' @title Rotated Copulas Created from an Existing Copula and a Mask of Logicals
##' @param copula The 'base' copula
##' @param flip A vector of logicals; if element i is TRUE,
##'        the copula is "rotated" wrt the axis x_i = 0.5;
##'        the default value is all TRUE which gives the survival copula
##' @return a new "rotCopula" object; see above
##' @author Ivan Kojadinovic
rotCopula <- function(copula, flip = rep(TRUE, copula@dimension)) {
    new("rotCopula",
        dimension = copula@dimension,
        parameters = copula@parameters,
        param.names = copula@param.names,
        param.lowbnd = copula@param.lowbnd,
        param.upbnd = copula@param.upbnd,
        copula = copula,
        flip = flip,
        fullname = paste("Rotated copula based on:", copula@fullname))
}

## Internal. swicth u[,i] to 1 - u[,i] according to flip
apply.flip <- function(u, flip) {
    if (!is.matrix(u)) u <- matrix(u, ncol = length(flip))
    u[,flip] <- 1 - u[,flip]
    u
}

## pCopula
pRotCopula <- function(u, copula) {
    copula@copula@parameters <- copula@parameters
    apply(apply.flip(u, copula@flip), 1, # TODO: vectorize prob ?
          function(x) prob(copula@copula,
                           l = pmin(x, copula@flip),
                           u = pmax(x, copula@flip)))
}

## dCopula
dRotCopula <- function(u, copula, log = FALSE, ...) {
    copula@copula@parameters <- copula@parameters
    dCopula(apply.flip(u, copula@flip), copula@copula, log = log, ...)
}

## rCopula
rRotCopula <- function(n, copula) {
    copula@copula@parameters <- copula@parameters
    apply.flip(rCopula(n, copula@copula), copula@flip)
}

## rho
rhoRotCopula <- function(copula) {
    ## if (copula@dimension > 2L)
    ##     warning("considering only the first bivariate margin of the copula")
    copula@copula@parameters <- copula@parameters
    (-1)^sum(copula@flip[1:2]) * rho(copula@copula)
}

## iRho
iRhoRotCopula <- function(copula, rho) {
    ## if (copula@dimension > 2L)
    ##     warning("considering only the first bivariate margin of the copula")
    iRho(copula@copula, (-1)^sum(copula@flip[1:2]) * rho)
}

## tau
tauRotCopula <- function(copula) {
    ## if (copula@dimension > 2L)
    ##     warning("considering only the first bivariate margin of the copula")
    copula@copula@parameters <- copula@parameters
    (-1)^sum(copula@flip[1:2]) * tau(copula@copula)
}

## iRho
iTauRotCopula <- function(copula, tau) {
    ## if (copula@dimension > 2L)
    ##     warning("considering only the first bivariate margin of the copula")
    iTau(copula@copula, (-1)^sum(copula@flip[1:2]) * tau)
}

## lambda
lambdaRotCopula <- function(copula) {
    ## if (copula@dimension > 2L)
    ##     warning("considering only the first bivariate margin of the copula")
    copula@copula@parameters <- copula@parameters
    sm <- sum(copula@flip[1:2])
    if (sm == 1) {
        warning("lambda() method for 'rotCopula' not implemented yet")
        c(lower= NA, upper= NA)
    }
    else {
        ti <- lambda(copula@copula)
        names(ti) <- NULL
        if (sm == 0)
            c(lower = ti[1], upper = ti[2])
        else
            c(lower = ti[2], upper = ti[1])
    }
}

## dRho
dRhoRotCopula <- function(copula) {
    ## if (copula@dimension > 2)
    ##     warning("considering only the first bivariate margin of the copula")
    copula@copula@parameters <- copula@parameters
    (-1)^sum(copula@flip[1:2]) * dRho(copula@copula)
}

## dTau
dTauRotCopula <- function(copula) {
    ## if (copula@dimension > 2)
    ##     warning("considering only the first bivariate margin of the copula")
    copula@copula@parameters <- copula@parameters
    (-1)^sum(copula@flip[1:2]) * dTau(copula@copula)
}

## dlogcdu
dlogcduRotCopula <- function(copula, u, ...) {
    copula@copula@parameters <- copula@parameters
    dlogcdu(copula@copula, apply.flip(u, copula@flip)) # TODO: check
}

## dCdu
dCduRotCopula <- function(copula, u, ...) {
    copula@copula@parameters <- copula@parameters
    dCdu(copula@copula, apply.flip(u, copula@flip)) # TODO: check
}

## dlogcdtheta
dlogcdthetaRotCopula <- function(copula, u, ...) {
    copula@copula@parameters <- copula@parameters
    dlogcdtheta(copula@copula, apply.flip(u, copula@flip)) # TODO: check
}

## dCdtheta
dCdthetaRotCopula <- function(copula, u, ...) {
    copula@copula@parameters <- copula@parameters
    dCdtheta(copula@copula, apply.flip(u, copula@flip)) # TODO: check
}

setMethod("pCopula", signature("numeric", "rotCopula"),pRotCopula)
setMethod("pCopula", signature("matrix", "rotCopula"), pRotCopula)

setMethod("dCopula", signature("numeric", "rotCopula"), dRotCopula)
setMethod("dCopula", signature("matrix", "rotCopula"), dRotCopula)

setMethod("rCopula", signature("numeric", "rotCopula"), rRotCopula)

setMethod("rho", signature("rotCopula"), rhoRotCopula)
setMethod("tau", signature("rotCopula"), tauRotCopula)

setMethod("iRho", signature("rotCopula"), iRhoRotCopula)
setMethod("iTau", signature("rotCopula"), iTauRotCopula)

setMethod("lambda", signature("rotCopula"), lambdaRotCopula)

setMethod("dRho", signature("rotCopula"), dRhoRotCopula)
setMethod("dTau", signature("rotCopula"), dTauRotCopula)

setMethod("dlogcdu", signature("rotCopula"), dlogcduRotCopula)
setMethod("dCdu", signature("rotCopula"), dCduRotCopula)
setMethod("dlogcdtheta", signature("rotCopula"), dlogcdthetaRotCopula)
setMethod("dCdtheta", signature("rotCopula"), dCdthetaRotCopula)
