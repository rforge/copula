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

##' @title A scatter plot matrix with nice variable names
##' @param data numeric matrix or as.matrix(.)able
##' @param varnames variable names, typically unspecified
##' @param Vname character string to become "root variable name"
##' @param col.mat matrix of colors
##' @param bg.col.mat matrix of background colors
##' @param ... further arguments to splom()
##' @return a splom() object
##' @author Martin Maechler
splom2 <- function(data, varnames=NULL, Vname="U", xlab = "",
                   col.mat = NULL, bg.col.mat = NULL, ...)
{
    stopifnot(require(lattice),
	      is.numeric(data <- as.matrix(data)),
	      (d <- ncol(data)) >= 1)
    if(is.null(varnames)) {
	varnames <- do.call(expression,
			    lapply(1:d, function(i)
				   substitute(italic(A[I]), list(A = as.name(Vname), I=0+i))))
    }
    n <- nrow(data)
    if(is.null(col.mat))
        col.mat <- matrix(trellis.par.get("plot.symbol")$col, n,d)
    if(is.null(bg.col.mat))
        bg.col.mat <- matrix(trellis.par.get("background")$col, n,d)
    ## From Deepayan Sarkar, working around missing feature
    ##		(which should be in next release) of lattice
    my.diag.panel <- function(x, varname, ...)
        diag.panel.splom(x, varname=parse(text=varname), ...)
    ## splom
    splom(~data[,1:d], varnames=varnames, diag.panel=my.diag.panel, xlab="",
          panel = function(x, y, i, j, ...) {
              panel.fill(bg.col.mat[i,j])
              panel.splom(x, y, col=col.mat[i,j])
          }, ...)
}


