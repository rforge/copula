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

####  'interval'	 class utilities
####  =========================== these are small and simple
###   use require(package= "Intervals")	 if you want serious interval "work"
interval <- function(ch) {
    ## Purpose: "interval" object constructor from string  "[ a, b)", ...
    ## Author: Martin Maechler, Date: 16 Nov 2009
    stopifnot(is.character(ch), length(ch) == 1L)
    sp <- strsplit(ch, ", *")[[1]]
    if(length(sp) != 2L) stop("'ch' must contain exactly one \",\"")
    L <- gsub(" +", "", sp[1]); bL <- substr(L, 1,1)
    if(!any(iL <- bL == c("(","[","]")))
	stop("interval specification must start with \"[\",  \"(\"  or \"]\"")
    R <- gsub(" +", "", sp[2]); nR <- nchar(R); bR <- substr(R, nR,nR)
    if(!any(iR <- bR == c(")","]","[")))
	stop("interval specification must end with \")\",  \"]\"  or \"[\"")
    new("interval", as.numeric(c(substring(L, 2), substr(R, 1, nR-1))),
	open = c(which(iL) != 2, which(iR) != 2))
}

setMethod("format", "interval",
	  function(x, trim = TRUE, ...) {
    r <- format(x@.Data, trim=trim, ...)
    paste(if(x@open[1]) "(" else "[", r[1],", ", r[2],
	  if(x@open[2]) ")" else "]", sep="")
})

setMethod("show", "interval",
	  function(object) cat("'interval' object  ", format(object), "\n",
	sep=''))

##' Summary group method: range(), min(), max(), [sum(), prod(), any(), all()] :
setMethod("Summary", signature(x = "interval", na.rm = "ANY"),
	  function(x, ..., na.rm) callGeneric(x@.Data, ..., na.rm=na.rm))
