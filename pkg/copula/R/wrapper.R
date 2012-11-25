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


### Wrappers and auxiliary functions for dealing with elliptical (Gauss, t_nu)
### and Archimedean copulas

##' Determine the copula class for a given copula object
##'
##' @title Copula class for the given copula object
##' @param cop copula object
##' @return "ellipCopula" or "outer_nacopula" depending on the given copula object
##' @author Marius Hofert
copClass <- function(cop)
{
    cls <- class(cop)
    if(is(cop, "copula") && (cls=="normalCopula" || cls=="tCopula")) "ellipCopula" # note: there could be other "copula" objects which are not elliptical
    else if(cls=="outer_nacopula") "outer_nacopula" # can be Archimedean or nested Archimedean
    else stop("not yet supported copula object")
}

##' Determine the copula family for a given copula object
##'
##' @title Copula family for the given copula object
##' @param cop copula object (either elliptical or (nested) Archimedean)
##' @return family string
##' @author Marius Hofert
copFamily <- function(cop)
{
    cls <- class(cop)
    if(is(cop, "copula")){
        if(cls=="normalCopula") "normal"
        else if(cls=="tCopula") "t"
        else stop("unsupported copula family")
    } else if(cls=="outer_nacopula"){
        cop@copula@name # could be nested or not
    } else stop("not yet supported copula object")
}

##' Determine the copula family for a given copula object
##'
##' @title Copula family for the given copula object
##' @param cop copula object (either elliptical or (nested) Archimedean)
##' @return family string
##' @author Marius Hofert
copFamilyClass <- function(family)
{
    if(family=="normal" || family=="t") "ellipCopula"
    else if(family %in% c_longNames ||
            family %in% paste0("opower:", c_longNames)) "outer_nacopula" # note: opower not really supported yet
    else stop("family ", family, " not yet supported")
}

##' Creating elliptical and Archimedean copula objects
##'
##' @title Creating elliptical and Archimedean copula objects
##' @param family specified family (for the elliptical cops, it's "normal" or "t")
##' @param theta parameter (a number or vector)
##' @param d dimension (can be omitted if nacList is given)
##' @param ... additional arguments:
##'        either 'dispstr', 'df', 'df.fixed' for the elliptical copulas
##'           or  'nacList' for the nested Archimedean {use 'theta' above for Archimedean}
##' @return a copula object
##' @author Marius Hofert and Martin Maechler
copCreate <- function(family, theta, d, ...)
{
    switch(copFamilyClass(family),
           "ellipCopula"={
               stopifnot(is.numeric(theta))
               ellipCopula(family, param=theta, dim=d, ...)
               ## Note: ... args are:
               ## dispstr="ex" (exchangeable), "ar1" (AR1), "toep" (Toeplitz), "un"
               ##         (unstructured)
               ## df (degrees of freedom)
               ## df.fixed=FALSE (TRUE means that d.o.f. will be considered as fixed
               ##          and won't be estimated)
           },
           "outer_nacopula"={
               L <- if(hasArg(nacList)){ # nested Archimedean case => we don't need theta
                   list(...)$nacList
               } else {
                   ## Archimedean case => one parameter
                   if(!is.numeric(theta) || length(theta)!=1)
                       stop("must specify 'nacList' (list)  or 'theta' (numeric)")
                   list(theta, 1:d)
               }
               stopifnot(is.list(L))
               onacopulaL(family, L)
               ## Note: ... with an argument "nacList" can be used to construct
               ## a nested Archimedean copula
           },
           stop("family ", family, " not yet supported"))
}

##' @title Construct a symmetric matrix with 1s on the diagonal from the given
##'        parameter vector
##' @param param parameter vector
##' @param d number of columns (or rows) of the output matrix
##' @return a symmetric matrix with 1s on the diagonal and the values of param
##'         filled column-wise below the diagonal (= row-wise above the diagonal)
##' @author Marius Hofert
p2P <- function(param, d){
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
##' Note: - This is used "by foot" at several points in the package.
##'       - A nice check is
##'         param <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
##'         stopifnot(P2p(p2P(param, 4))==param)
P2p <- function(P) P[lower.tri(P)]

