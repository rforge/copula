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
        fullname = paste("Rotated copula based on: [", copula@fullname, "]"))
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
    copula@copula@parameters <- copula@parameters
    (-1)^sum(copula@flip[1:2]) * rho(copula@copula)
}

## iRho
iRhoRotCopula <- function(copula, rho) {
    iRho(copula@copula, (-1)^sum(copula@flip[1:2]) * rho)
}

## tau
tauRotCopula <- function(copula) {
    copula@copula@parameters <- copula@parameters
    (-1)^sum(copula@flip[1:2]) * tau(copula@copula)
}

## iRho
iTauRotCopula <- function(copula, tau) {
    iTau(copula@copula, (-1)^sum(copula@flip[1:2]) * tau)
}

## lambda
lambdaRotCopula <- function(copula) {
    ## if (copula@dimension > 2L)
    ##     warning("considering only the first bivariate margin of the copula")
    copula@copula@parameters <- copula@parameters
    sm <- sum(copula@flip[1:2])
    if (sm == 1) {
        warning("lambda() method for copula class 'rotCopula' not implemented yet")
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

## fitCopula
fitCopulaRotCopula  <- function(copula, data, ...) {
    if(!is.matrix(data)) {
        warning("coercing 'data' to a matrix.")
        data <- as.matrix(data); stopifnot(is.matrix(data))
    }
    fit <- fitCopulaCopula(copula@copula, data = apply.flip(data, copula@flip), ...)
    copula@copula <- fit@copula
    copula@parameters <- fit@estimate
    fit@copula <- copula
    fit
}

## gofCopula
gofCopulaRotCopula  <- function(copula, x, ...) {
    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    x[,copula@flip] <- -x[,copula@flip]
    gofCopulaCopula(copula@copula, x = x, ...)
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

setMethod("fitCopula", signature("rotCopula"), fitCopulaRotCopula)
setMethod("gofCopula", signature("rotCopula"), gofCopulaRotCopula)
