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


mvdc <- function(copula, margins, paramMargins, marginsIdentical = FALSE) {
  if (marginsIdentical) {
    if (length(margins) == 1) margins <- rep(margins, copula@dimension)
    if (length(paramMargins) == 1) paramMargins <- rep(paramMargins, copula@dimension)
  }
  new("mvdc", copula = copula, margins = margins, paramMargins = paramMargins,
       marginsIdentical = marginsIdentical)
}

##' @title Parameter names of the margins of an "mvdc" object
##' @param mv 
##' @return character vector of "the correct" length
##' @author Martin Maechler
margpnames <- function(mv) {
    nMar <- vapply(mv@paramMargins, length, 1L)
    p <- mv@copula@dimension
    pnms <- unlist(lapply(mv@paramMargins, names)) # maybe NULL
    if (sum(nMar) == 0) character()
    else if(mv@marginsIdentical)
	paste(paste("m", pnms[1], sep="."))
    else
	paste(paste0("m", rep.int(1:p, nMar)), pnms, sep=".")
}

## Function asCall was kindly supplied by
## Martin Maechler <maechler@stat.math.ethz.ch>,
## motivated by an application of nor1mix and copula
## from Lei Liu <liulei@virginia.edu>.
## They fixes the function getExpr in the old
## version, which assumed that the parameters to
## [rdpq]<distrib> were vectors.

asCall <- function(fun, param)
{
    cc <-
	if (length(param) == 0)
	    quote(FUN(x))
	else if(is.list(param)) {
	    as.call(c(quote(FUN), c(quote(x), as.expression(param))))
	} else { ## assume that [dpq]<distrib>(x, param) will work
	    as.call(c(quote(FUN), c(quote(x), substitute(param))))
	}
    cc[[1]] <- as.name(fun)
    cc
}

dMvdc <- function(x, mvdc, log=FALSE) {
  dim <- mvdc@copula@dimension
  densmarg <- if(log) 0 else 1
  if (is.vector(x)) x <- matrix(x, nrow = 1)
  u <- x
  for (i in 1:dim) {
    cdf.expr <- asCall(paste0("p", mvdc@margins[i]), mvdc@paramMargins[[i]])
    pdf.expr <- asCall(paste0("d", mvdc@margins[i]), mvdc@paramMargins[[i]])
    u[,i] <- eval(cdf.expr, list(x = x[,i]))
    densmarg <-
	if(log)
	    ## FIXME: user should be able to give density which has a log argument
	    densmarg + log(eval(pdf.expr, list(x = x[,i])))
	else
	    densmarg * eval(pdf.expr, list(x = x[,i]))
  }
  if(log)
      dCopula(u, mvdc@copula, log=TRUE) + densmarg
  else
      dCopula(u, mvdc@copula) * densmarg
}

pMvdc <- function(x, mvdc) {
  dim <- mvdc@copula@dimension
  if (is.vector(x)) x <- matrix(x, nrow = 1)
  u <- x
  for (i in 1:dim) {
    cdf.expr <- asCall(paste0("p", mvdc@margins[i]), mvdc@paramMargins[[i]])
    u[,i] <- eval(cdf.expr, list(x = x[,i]))
  }
  pCopula(u, mvdc@copula)
}

rMvdc <- function(n, mvdc) {
  dim <- mvdc@copula@dimension
  u <- rCopula(n, mvdc@copula)
  x <- u
  for (i in 1:dim) {
    qdf.expr <- asCall(paste0("q", mvdc@margins[i]), mvdc@paramMargins[[i]])
    x[,i] <- eval(qdf.expr, list(x = u[,i]))
  }
  x
}

dmvdc <- function(mvdc, x, log=FALSE) { .Deprecated("dMvdc"); dMvdc(x, mvdc, log) }
pmvdc <- function(mvdc, x) { .Deprecated("pMvdc"); pMvdc(x, mvdc) }
rmvdc <- function(mvdc, n) { .Deprecated("rMvdc"); rMvdc(n, mvdc) }

