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


archmCopula <- function(family, param = NA_real_, dim = 2L, ...) {
    fams <- sub("Copula$", '', names(getClass("archmCopula")@subclasses))
    fam <- pmatch(family <- tolower(family), fams)
    if(is.na(fam))
        stop("Valid family names are ", paste(dQuote(fams), collapse=", "))
    dim <- as.integer(dim)
    if(family == "amh" && dim != 2L)
	stop("'amh' is not yet available for dim > 2")
    switch(fam,
	   claytonCopula(param, dim = dim),
	   frankCopula	(param, dim = dim),
	   amhCopula	(param, dim = dim),
	   gumbelCopula (param, dim = dim)
	   )
}


kendallsTauArchmCopula <- function(copula) {
  integrand <- function(x) genFun(copula, x) / genFunDer1(copula, x)
  1 + 4 * integrate(integrand, 0, 1)$value
}

setMethod("kendallsTau", signature("archmCopula"), kendallsTauArchmCopula)
