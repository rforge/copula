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


##' Empirical copula of U at u
##'
##' @title Empirical copula of U at u
##' @param u (m, d) matrix of evaluation points
##' @param U (n, d) matrix of pseudo-data based on which the empirical copula
##'        is computed (if not pseudo-data already, use do.pobs=TRUE)
##' @param do.pobs logical indicating whether pobs is applied to U
##' @param offset scaling factor sum()/(n+offset) when computing the empirical
##'        copula
##' @param method method string
##' @return empirical copula of U at u
##' @author Ivan Kojadinovic and Marius Hofert
##' Note: See the .Rd for a nice graphical check with the Kendall function
Cn <- function(u, U, do.pobs=TRUE, offset=0, method=c("C", "R"))
{
    if(!is.matrix(u)) u <- cbind(u, deparse.level=0L)
    if(!is.matrix(U)) U <- cbind(U, deparse.level=0L)
    stopifnot((d <- ncol(U))==ncol(u),
              0 <= u, u <= 1, 0 <= U, U <= 1)
    if(do.pobs) U <- pobs(U)
    n <- nrow(U)

    ## d = 1
    if(d==1) return( vapply(u, function(u.) sum(U<=u.)/(n+offset), numeric(1)) )

    ## d > 1
    method <- match.arg(method)
    switch(method,
           "C"={
               m <- nrow(u)
               .C(RmultCn, # see empcop.c
                  as.double(U),
                  as.integer(n),
                  as.integer(d),
                  as.double(u),
                  as.integer(m),
                  ec=double(m),
                  as.double(offset))$ec
           },
           "R"={
               apply(u, 1, function(u.) sum(colSums(t(U)<=u.)==d)/(n+offset) )
           },
           stop("wrong method"))
}
