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


##' @title Conditional copula
##' @param u data matrix in [0,1]^(n, d) ((pseudo-/copula-)observations if
##'        inverse==TRUE and U[0,1] observations if inverse==FALSE)
##' @param cop object of class Copula
##' @param j.ind single index in {2,..,d} for which C(u_j | u_1,..,u_{j-1}) is
##'        computed and returned or NULL for which C(u_j | u_1,..,u_{j-1}) is
##'        computed for all j in {2,...,d} (and in which case the non-transformed
##'        first column is also returned)
##' @param n.MC parameter n.MC for evaluating the derivatives via Monte Carlo;
##'        if 0, the available (theoretical) formula is used.
##' @param inverse logical indicating whether the inverse of rtrafo is computed
##'        (this is known as 'conditional distribution method' for sampling)
##' @param log logical indicating whether the log-transform is computed
##' @return An (n,d) matrix U of supposedly U[0,1]^d realizations (if inverse==FALSE)
##'         or copula samples (if inverse=TRUE) (if j.ind==NULL) or
##'         C(u_j | u_1,..,u_{j-1}) (if j.ind in {2,..,d}) [or the log of
##'         the result if log=TRUE]
##' @author Marius Hofert and Martin Maechler
##' @note See ./gofTrafos.R
cCopula <-  function(u, copula, j.ind=ncol(u), n.MC=0, inverse=FALSE, log=FALSE)
{
    stopifnot(is(copula, "Copula"))
    drop(rtrafo(u, cop=copula, j.ind=j.ind, n.MC=n.MC, inverse=inverse, log=log))
}

