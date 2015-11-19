## Copyright (C) 2012 Marius Hofert and Martin Maechler
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


### Tools around matrices

##' @title Construct a symmetric matrix with 1s on the diagonal from the given
##'        parameter vector
##' @param param parameter vector
##' @param d number of columns (or rows) of the output matrix
##' @return a symmetric matrix with 1s on the diagonal and the values of param
##'         filled column-wise below the diagonal (= row-wise above the diagonal)
##' @author Marius Hofert
p2P <- function(param, d)
{
    P <- diag(1, nrow=d)
    P[lower.tri(P)] <- param
    P <- P+t(P)
    diag(P) <- rep.int(1, d)
    P
}

##' @title Extract the vector of column-wise below-diagonal entries from a matrix
##' @param P matrix (typically a symmetric matrix as used for elliptical copulas)
##' @return the vector of column-wise below-diagonal entries of P (they are equal
##'         to the row-wise above-diagonal entries in case of a symmetric matrix)
##' @author Marius Hofert
##' Note: This is used "by foot" at several points in the package.
P2p <- function(P) P[lower.tri(P)]

##' @title Construct matrix Sigma from a given elliptical copula
##' @param copula copula
##' @return (d, d) matrix Sigma containing the parameter vector rho
##' @author Marius Hofert
getSigma <- function(copula)
{
    stopifnot(is(copula, "ellipCopula"))
    d <- copula@dimension
    rho <- copula@getRho(copula)
    switch(copula@dispstr,
	   "ex" = {
	       Sigma <- matrix(rho[1], nrow=d, ncol=d)
	       diag(Sigma) <- rep(1, d)
	       Sigma
	   },
	   "ar1" = {
	       rho^abs(outer(1:d, 1:d, FUN="-"))
	   },
	   "un" = {
	       p2P(rho, d)
	   },
	   "toep" = {
	       rho <- c(rho, 1)
	       ind <- outer(1:d, 1:d, FUN=function(i, j) abs(i-j))
	       diag(ind) <- length(rho)
	       matrix(rho[ind], nrow=d, ncol=d)
	   },
	   stop("invalid 'dispstr'"))
}

##' @title Find the pairs with smallest (or largest) n values in a symmetric
##'        matrix
##' @param x A symmetric matrix
##' @param n Number of extreme values to be returned
##' @param decreasing logical indicating whether the numbers are sorted in
##'        decreasing order
##' @return A (n, 3)-matrix with the smallest (or, if decreasing=TRUE, largest)
##'         n values in the symmetric matrix x (3rd col) and the corresponding
##'         indices (1st and 2nd col)
##' @author Marius Hofert
extreme_pairs <- function(x, n=6, decreasing=FALSE)
{
    if(!is.matrix(x)) x <- rbind(x, deparse.level=0L)
    d <- ncol(x)
    stopifnot(n>=1, is.logical(decreasing), d>=2, nrow(x)==d)
    ind <- expand.grid(1:d, 1:d)
    names(ind) <- c("row", "col")
    ind <- ind[ind[,"row"]>ind[,"col"],] # pick out indices of lower triangular matrix
    x.vec <- P2p(x) # values of lower triangular matrix as a vector
    val.ind <- cbind(ind, value=x.vec)
    val.ind.sorted <- val.ind[order(val.ind[,"value"], decreasing=decreasing),] # sort according to values
    rownames(val.ind.sorted) <- NULL
    val.ind.sorted[1:min(n, d*(d-1)/2),]
}
