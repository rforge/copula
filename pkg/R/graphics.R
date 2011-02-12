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
                   col.mat=matrix(rep("black", ncol(data)^2), ncol=ncol(data)), 
                   bg.col.mat=matrix(rep(NA, ncol(data)^2), ncol=ncol(data)), ...)
{
    stopifnot(require(lattice),
	      is.numeric(data <- as.matrix(data)),
	      (d <- ncol(data)) >= 1)
    if(is.null(varnames)) {
	varnames <- do.call(expression,
			    lapply(1:d, function(i)
				   substitute(A[I], list(A = as.name(Vname), I=0+i))))
    }
    ## From Deepayan Sarkar, working around missing feature
    ##		(which should be in next release) of lattice
    my.diag.panel <- function(x, varname, ...)
        diag.panel.splom(x, varname=parse(text=varname), ...)
    ## splom
    splom(~data[,1:d], varnames=varnames, diag.panel=my.diag.panel, xlab="",
          panel=function(x, y, i, j, ...){
              panel.fill(bg.col.mat[i,j])
              panel.splom(x, y, col=col.mat[i,j])
          }, ...)
}

##' Plots a scatter plot matrix of the provided data
##'
##' @title Adjusted scatter plot matrix
##' @param data data matrix
##' @param device graphic device to be used - as in trellis.device()
##' @param color logical indicating if the plot is colored (as in trellis.device)
##' @param outfilename name of the output file (without file ending)
##' @param varnames variable names to be printed on the diagonal
##' @param ... additional arguments passed to the splom call
##' @return the lattice / grid plot object, invisibly
##' @author Marius Hofert, Martin Maechler
splomFOO <- function(data, device = getOption("device"),
                     color = !(dev.name == "postscript"),
                     varnames = NULL, Vname = "U", outfilename = "splom2", ...)
{
    stopifnot(require(lattice),
	      is.numeric(data <- as.matrix(data)),
	      (d <- ncol(data)) >= 1, # numeric matrix
	      is.character(outfilename))

    dev.name <-
        if (is.character(device)) device else deparse(substitute(device))
    if(is.null(varnames)) {
	varnames <- do.call(expression,
			    lapply(1:d, function(i)
				   substitute(A[I], list(A = as.name(Vname), I=0+i))))
    }
    isdeviceFile <- dev.name %in% c("pdf", "postscript", "png")

    ## FIXME: isn't that solved now?
    ## AAARGH: splom() will *NOT* work with expression varnames
    ## but the simple pairs() actually does:
    ## pairs(data, varNames, gap=0)   # ok

    if(isdeviceFile) {
        file <- paste(outfilename,
                      switch(device,
                             pdf = "pdf",
                             postscript = "ps",
                             png = "png"),
                      sep = ".")
        trellis.device(device = device, color = color, file = file)
    }
    else
        trellis.device(device = device, color = color)

    print(G <- splom(~data[,1:d], varnames = varnames, ...))

    if(isdeviceFile) {
        cat("closing trellis.device",.Device, "\n")
        dev.off()
    }
    invisible(G)
}

