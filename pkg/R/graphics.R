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

##' Plots a scatterplot matrix of the provided data
##' @param data data matrix
##' @param device
##' @param color bool if the plot is colored
##' @param outfilename name of the outfile (without file ending)
##' @param varnames variable names to be printed on the diagonal
##' @param ... additional arguments passed to the splom call
##' @return none
##' @author Marius Hofert, Martin Maechler
require(lattice)
spm <- function(data, dev = c("pdf", "postscript", "png"), 
                color = !(dev.name == "postscript"), outfilename = "plot", ...){
    data <- as.matrix(data)
    stopifnot(is.matrix(data), # check data coercion  
              all(dev %in% c("pdf", "postscript", "png")), # check device 
              all(color %in% c(TRUE,FALSE)), # check color
              is.character(outfilename) # check outfile name
              ) 
    
    d <- dim(data)[2] # get the dimension of the data
    ## if(varnames argument to splom is not provided){ # build default vector of varnames
    vec <- "U[1]"
    for(i in 2:d) vec <- c(vec,paste("U[",i,"]",sep=""))
    varnames <- expression(vec)
    ## }
    outfile <- getwd() # build outfile path
    outfile <- paste(outfile,outfilename,sep="")
    if(dev == "pdf"){
	outfile <- paste(outfile,".pdf",sep="")
    } else if(dev == "postscript"){
	outfile <- paste(outfile,".ps",sep="")
    } else {
	outfile <- paste(outfile,".png",sep="")
    }
    ## call splom
    trellis.device(device = dev, color = color,file = outfile, 
                   useKerning = FALSE)
    print(splom(~data[,1:d], ...))
    dev.off()
}

