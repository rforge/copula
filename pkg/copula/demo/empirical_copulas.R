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
##
## Implementing various (smoothed) empirical copulas in R


### Setup ######################################################################

library(copula)


### Auxiliary functions ########################################################

##' @title Compute log(x_1 + .. + x_n) for given log(x_1),..,log(x_n)
##' @param lx n-vector containing log(x_1),..,log(x_n)
##' @author Marius Hofert
lsum <- function(lx) {
    l.off <- max(lx)
    l.off + log(sum(exp(lx - l.off)))
}

##' @title Evaluate Empirical Copulas
##' @param u evaluation points (in [0,1]^d)
##' @param U points (in [0,1]^d) based on which the empirical copula is built
##' @param smoothing character string indicating the smoothing method
##' @param offset shift when computing the outer sum (over i)
##' @param log logical indicating whether the log-density is to be returned
##'        (only applies to smoothing = "dbeta")
##' @param ... additional arguments passed to the underlying rank()
##' @return empirical copula values
##' @author Marius Hofert
##' @note See Hofert et al. (2018, "Elements of copula modeling with R")
empirical_copula <- function(u, U, smoothing = c("none", "pbeta", "dbeta", "checkerboard"),
                             offset = 0, log = FALSE, ...)
{
    stopifnot(0 <= u, u <= 1, 0 <= U, U <= 1)
    if(!is.matrix(u)) u <- rbind(u)
    if(!is.matrix(U)) U <- rbind(U)
    m <- nrow(u) # number of evaluation points
    n <- nrow(U) # number of points based on which the empirical copula is computed
    d <- ncol(U) # dimension
    stopifnot(ncol(u) == d)
    R <- apply(U, 2, rank, ...) # (n, d)-matrix of ranks
    switch(match.arg(smoothing),
    "none" = {
        R. <- t(R) / (n + 1) # (d, n)-matrix
        vapply(seq_len(m), function(k) # iterate over rows k of u
            sum(colSums(R. <= u[k,]) == d) / (n + offset), NA_real_)
    },
    "pbeta" = {
        ## Note: pbeta(q, shape1, shape2) is vectorized in the following sense:
        ##       1) pbeta(c(0.8, 0.6), shape1 = 1, shape2 = 1) => 2-vector (as expected)
        ##       2) pbeta(0.8, shape1 = 1:4, shape2 = 1:4) => 4-vector (as expected)
        ##       3) pbeta(c(0.8, 0.6), shape1 = 1:2, shape2 = 1:2) => 2-vector (as expected)
        ##       4) pbeta(c(0.8, 0.6), shape1 = 1:4, shape2 = 1:4) => This is
        ##          equal to the recycled pbeta(c(0.8, 0.6, 0.8, 0.6), shape1 = 1:4, shape2 = 1:4)
        vapply(seq_len(m), function(k) { # iterate over rows k of u
                sum( # sum() over i
                    vapply(seq_len(n), function(i)
                        prod( pbeta(u[k,], shape1 = R[i,], shape2 = n + 1 - R[i,]) ), # prod() over j
                        NA_real_)) / (n + offset)
        }, NA_real_)
    },
    "dbeta" = {
        if(log) {
            vapply(seq_len(m), function(k) { # iterate over rows k of u
                lsum( # lsum() over i
                    vapply(seq_len(n), function(i) {
                        ## k and i are fixed now
                        lx.k.i <- sum( dbeta(u[k,], shape1 = R[i,], shape2 = n + 1 - R[i,], log = TRUE) ) # log(prod()) = sum(log()) over j for fixed k and i
                    },
                    NA_real_)) - log(n + offset)
            }, NA_real_)
        } else { # as for 'pbeta', just with dbeta()
            vapply(seq_len(m), function(k) { # iterate over rows k of u
                sum( # sum() over i
                    vapply(seq_len(n), function(i)
                        prod( dbeta(u[k,], shape1 = R[i,], shape2 = n + 1 - R[i,]) ), # prod() over j
                        NA_real_)) / (n + offset)
            }, NA_real_)
        }
    },
    "checkerboard" = {
        R. <- t(R) # (d, n)-matrix
        vapply(seq_len(m), function(k) # iterate over rows k of u
            sum(apply(pmin(pmax(n * u[k,] - R. + 1, 0), 1), 2, prod)) / (n + offset),
            NA_real_) # pmin(...) = (d, n)-matrix
    },
    stop("Wrong 'smoothing'"))
}


### Main #######################################################################

## Generate copula data based on which to build empirical copula
n <- 1000 # sample size
d <- 2 # dimension
set.seed(271)
U <- rCopula(n, copula = gumbelCopula(iTau(gumbelCopula(), tau = 0.5), dim = d))

## Evaluation points
n.grid <- 26
sq <- seq(0, 1, length.out = n.grid)
u <- as.matrix(expand.grid("u[1]" = sq, "u[2]" = sq, KEEP.OUT.ATTRS = FALSE))
## ... for the density of the empirical beta copula
delta <- 0.01
sq. <- seq(delta, 1-delta, length.out = n.grid)
u. <- as.matrix(expand.grid("u[1]" = sq., "u[2]" = sq., KEEP.OUT.ATTRS = FALSE))

## Evaluate empirical copulas
emp.cop.none   <- empirical_copula(u,  U = U)
emp.cop.pbeta  <- empirical_copula(u,  U = U, smoothing = "pbeta")
emp.cop.dbeta  <- empirical_copula(u., U = U, smoothing = "dbeta")
lemp.cop.dbeta <- empirical_copula(u., U = U, smoothing = "dbeta", log = TRUE)
stopifnot(all.equal(lemp.cop.dbeta, log(emp.cop.dbeta))) # sanity check
emp.cop.chck   <- empirical_copula(u,  U = U, smoothing = "checkerboard")

## Plot and comparison with the functions of 'copula'
## Empirical copula
wireframe2(cbind(u, "Empirical~copula" = emp.cop.none))
wireframe2(empCopula(U), FUN = pCopula, draw.4.pCoplines = FALSE)
## Empirical beta copula
wireframe2(cbind(u, "Empirical~'beta'~copula" = emp.cop.pbeta))
wireframe2(empCopula(U, smoothing = "beta"), FUN = pCopula, draw.4.pCoplines = FALSE)
## Empirical beta copula density
wireframe2(cbind(u., "Empirical~'beta'~copula~density" = emp.cop.dbeta))
wireframe2(empCopula(U, smoothing = "beta"), FUN = dCopula, delta = delta, draw.4.pCoplines = FALSE)
## Empirical checkerboard copula
wireframe2(cbind(u, "Empirical~checkerboard~copula" = emp.cop.chck))
wireframe2(empCopula(U, smoothing = "checkerboard"), FUN = pCopula, draw.4.pCoplines = FALSE)
