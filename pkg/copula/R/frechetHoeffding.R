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


### Lower Frechet--Hoeffding bound W ###########################################

## Constructor
W <- function(dim = 2L)
{
    if(dim != 2) stop("W only available in the bivariate case")
    cdfExpr <- function(d) {
        uis <- c(paste0("u", 1:d, "-1"), "1")
        expr <- paste(uis, collapse = "+")
        expr <- paste0("max(", expr, ",0)") # TODO: only for a single 'u' okay? otherwise we would need 'pmax'
        parse(text = expr)
    }
    cdf <- cdfExpr(dim)
    new("W",
        dimension = dim,
        exprdist = c(cdf = cdf, # TODO: needed for what?
                     pdf = expression()), # TODO: no density => okay like this?
        parameters = double(0),
        param.names = character(0),
        param.lowbnd = double(0),
        param.upbnd = double(0),
        fullname = "<deprecated slot>")
}


### Methods ####################################################################

setMethod("rCopula", signature("numeric", "W"),
          function(n, copula) {
              U <- runif(n)
              cbind(U, 1-U)
          })
setMethod("pCopula", signature("matrix", "W"),
	  function(u, copula, log.p = FALSE) {
              d <- ncol(u)
              stopifnot(d == copula@dimension)
              res <- pmax(rowSums(u) - d + 1, 0)
	      if(log.p) log(res) else res
	  })
setMethod("dCopula", signature("matrix", "W"),
	  function(u, copula, log = FALSE, ...) {
	      stopifnot(ncol(u) == copula@dimension)
	      rep.int(if(log) -Inf else 0, nrow(u)) # TODO: okay? technically not a density, might be safer to return an error (?)
	  })

setMethod("tau", "W", function(copula, ...) -1)
setMethod("rho", "W", function(copula, ...) -1)
setMethod("lambda", "W", function(copula, ...)  c(lower = 0, upper = 0))

## GETR
setMethod("describeCop", c("W", "character"),
          function(x, kind = c("short", "very short", "long"), prefix = "", ...) {
    kind <- match.arg(kind)
    if(kind == "very short") # e.g. for show() which has more parts
        return(paste0(prefix, "W copula"))
    ## else
    d <- dim(x)
    ch <- paste0(prefix, "W copula, dim. d = ", d)
    switch(kind <- match.arg(kind),
           short = ch,
           long = ch,
           stop("invalid 'kind': ", kind))
})


### Upper Frechet--Hoeffding bound M ###########################################

## Constructor
M <- function(dim = 2L)
{
    cdfExpr <- function(d) {
        uis <- paste0("u", 1:d, collapse = ",")
        expr <- paste0("min(", uis, ")") # TODO: only for a single 'u' okay? otherwise we would need 'pmax'
        parse(text = expr)
    }
    cdf <- cdfExpr(dim)
    new("M",
        dimension = dim,
        exprdist = c(cdf = cdf, # TODO: needed for what?
                     pdf = expression()), # TODO: no density => okay like this?
        parameters = double(0),
        param.names = character(0),
        param.lowbnd = double(0),
        param.upbnd = double(0),
        fullname = "<deprecated slot>")
}


### Methods ####################################################################

setMethod("rCopula", signature("numeric", "M"),
          function(n, copula) matrix(runif(n), nrow = n, ncol = copula@dimension))
setMethod("pCopula", signature("matrix", "M"),
	  function(u, copula, log.p = FALSE) {
              stopifnot(ncol(u) == copula@dimension)
              res <- apply(u, 1, min)
	      if(log.p) log(res) else res
	  })
setMethod("dCopula", signature("matrix", "M"),
	  function(u, copula, log = FALSE, ...) {
	      stopifnot(ncol(u) == copula@dimension)
	      rep.int(if(log) -Inf else 0, nrow(u)) # TODO: okay? technically not a density, might be safer to return an error (?)
	  })

setMethod("tau", "M", function(copula, ...) 1)
setMethod("rho", "M", function(copula, ...) 1)
setMethod("lambda", "M", function(copula, ...)  c(lower = 1, upper = 1))

## GETR
setMethod("describeCop", c("M", "character"),
          function(x, kind = c("short", "very short", "long"), prefix = "", ...) {
    kind <- match.arg(kind)
    if(kind == "very short") # e.g. for show() which has more parts
        return(paste0(prefix, "M copula"))
    ## else
    d <- dim(x)
    ch <- paste0(prefix, "M copula, dim. d = ", d)
    switch(kind <- match.arg(kind),
           short = ch,
           long = ch,
           stop("invalid 'kind': ", kind))
})
