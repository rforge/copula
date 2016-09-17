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

##################################################################################
### Virtual class of asymmetric copulas
##################################################################################

setClass("asymCopula", contains = c("parCopula", "VIRTUAL"))

##################################################################################
### Virtual class of asymmetric copulas constructed from two d-dimensional copulas
##################################################################################

setClass("asym2Copula", contains = c("asymCopula", "VIRTUAL"),
	 slots = c(
             copula1 = "parCopula",
             copula2 = "parCopula",
             shapes = "numeric"
         ),
	 validity = function(object) {
    if(dim(object@copula1) != dim(object@copula2))
        "The argument copulas are not of the same dimension"
    else if(dim(object@copula1) != length(object@shapes))
        "The length of 'shapes' is not equal to the dimension of the copulas"
    else TRUE
})

##################################################################################
### Potentially asymmetric d-dimensional copulas of the form C(u^{1-a}) * D(u^a)
### This construction is known as Khoudraji's device
##################################################################################

setClass("khoudrajiCopula", contains = c("asym2Copula"))
## slots *and* validity currently from 'asym2Copula' !

##################################################################################
### Bivariate Khoudraji copulas
### Can be constructed from any bivariate copulas
##################################################################################

setClass("khoudrajiBivCopula", contains = "khoudrajiCopula",
	 validity = function(object)
    if(dim(object@copula1) != 2)
        "The argument copulas must be of dimension two"
    else TRUE
    )

################################################################################
### Khoudraji *explicit* copulas
### Explicit refers to the fact that the two d-dimensional copulas
### have explicit cdfs and pdfs
################################################################################

setClass("khoudrajiExplicitCopula", contains = "khoudrajiCopula",
	 slots = c(
             exprdist = "expression",
             derExprs1 = "expression",
             derExprs2 = "expression"
         )
         ## validity = function(object) {
         ##     ## TODO: check exprdist, derExprs[12]
         ## },
         )
##################################################################################
### Basic methods
##################################################################################

## dimension
setMethod("dim", signature("khoudrajiCopula"), function(x) dim(x@copula1))

## logical indicating which parameters are free
setMethod("isFree", signature("khoudrajiCopula"), function(copula)
    c(isFree(copula@copula1), isFree(copula@copula2), isFreeP(copula@shapes)))

## number of (free / all) parameters :
setMethod("nParam", signature("khoudrajiCopula"), function(copula, freeOnly=FALSE)
    nParam(copula@copula1, freeOnly=freeOnly) +
    nParam(copula@copula2, freeOnly=freeOnly) +
    (if(freeOnly) nFree else length)(copula@shapes))

## parameter names
setMethod("paramNames", signature("khoudrajiCopula"), function(x) {
    c(if (nParam(x@copula1, freeOnly=TRUE) > 0L) paste0("c1.", paramNames(x@copula1)) else NULL,
      if (nParam(x@copula2, freeOnly=TRUE) > 0L) paste0("c2.", paramNames(x@copula2)) else NULL,
      paste0("shape", 1:dim(x))[isFreeP(x@shapes)])
})

## get parameters
setMethod("getParam", signature("khoudrajiCopula"),
          function(copula, freeOnly = TRUE, attr = FALSE, named = attr) {
    par1 <- getParam(copula@copula1, freeOnly = freeOnly, attr = attr, named = named)
    par2 <- getParam(copula@copula2, freeOnly = freeOnly, attr = attr, named = named)
    fixed <- attr(copula@shapes, "fixed")
    d <- dim(copula)
    ns <- if (freeOnly) nFree(copula@shapes) else d
    sel <- if (!is.null(fixed) && freeOnly) !fixed else TRUE
    par <- if (named) c(par1, par2, setNames(copula@shapes[sel], paste0("shape", 1:d)[sel]))
           else c(par1, par2, copula@shapes[sel])
    if (!attr)
        par
    else
        structure(par,
                  param.lowbnd = c(attr(par1, "param.lowbnd"),
                                   attr(par2, "param.lowbnd"), rep(0, ns)),
                  param.upbnd  = c(attr(par1, "param.upbnd"),
                                   attr(par2, "param.upbnd"), rep(1, ns)))

})

## set free parameters
setMethod("freeParam<-", signature("khoudrajiCopula", "numeric"),
          function(copula, value) {
    n1 <- nParam(copula@copula1, freeOnly=TRUE)
    n2 <- nParam(copula@copula2, freeOnly=TRUE)
    ns <- nFree(copula@shapes)
    if (n1 + n2 + ns != length(value))
        stop("the length of 'value' is not equal to the number of free parameters")
    if (n1 > 0L) freeParam(copula@copula1) <- value[1:n1]
    if (n2 > 0L) freeParam(copula@copula2) <- value[n1 + 1:n2]
    if (ns > 0L) {
        fixed <- !isFreeP(copula@shapes)
        copula@shapes[!fixed] <- value[n1 + n2 + 1:ns]
    }
    copula
})

## set or modify "fixedness" of parameters
setMethod("fixedParam<-", signature("khoudrajiCopula", "logical"),
          function(copula, value) {
    stopifnot(length(value) %in% c(1L, nParam(copula)))
    if (identical(value, FALSE) || !any(value))
        copula
    if (anyNA(getParam(copula, freeOnly = FALSE)[value])) stop("Fixed parameters cannot be NA.")
    n1 <- nParam(copula@copula1)
    n2 <- nParam(copula@copula2)
    if (identical(value, FALSE) || !any(value)) {
         if (n1 > 0L) fixedParam(copula@copula1) <- FALSE
         if (n2 > 0L) fixedParam(copula@copula2) <- FALSE
         attr(copula@shapes, "fixed") <- NULL
    } else if (identical(value, TRUE) || all(value)) {
         if (n1 > 0L) fixedParam(copula@copula1) <- TRUE
         if (n2 > 0L) fixedParam(copula@copula2) <- TRUE
         attr(copula@shapes, "fixed") <- TRUE
    } else {
        if (n1 > 0L) fixedParam(copula@copula1) <- value[1:n1]
        if (n2 > 0L) fixedParam(copula@copula2) <- value[n1 + 1:n2]
        ns <- length(copula@shapes)
        attr(copula@shapes, "fixed") <-
            if (!any(value[n1 + n2 + 1:ns])) NULL
            else if (all(value[n1 + n2 + 1:ns])) TRUE else value[n1 + n2 + 1:ns]
    }
    copula
})

## describe copula
setMethod(describeCop, c("khoudrajiCopula", "character"), function(x, kind) {
    paste0("Khoudraji copula constructed from\n", describeCop(x@copula1, kind),
          "\nand\n", describeCop(x@copula2, kind))
})

##################################################################################
### Generators for Khoudraji copulas
##################################################################################

## C(u_1^{1-a_1}, u_2^{1-a_2}) * D(u_1^a_1, u_2^a_1) = C(g(u, 1-a)) * D(g(u, a))
KhoudFn <-
    list(
        g = function(u, a) u^a,
        ## inverse of g :
        ig = function(u, a) u^(1 / a),
        ## derivative wrt u :
        dgdu = function(u, a) a * u^(a - 1)
        )

##################################################################################
### Constructor of Khoudraji copulas
##################################################################################

##' Creates a khoudrajiBivCopula object, a khoudrajiExplicitCopula
##' or a khoudrajCopula object
##'
##' @title Creates a khoudrajiBivCopula object, a khoudrajiExplicitCopula
##' or a khoudrajCopula object
##' @param copula1 a copula
##' @param copula2 a copula
##' @param shapes a numeric of length dim(copula) with elements in [0,1]
##' @return a new khoudrajiBivCopula, khoudrajiExplicitCopula
##' or a khoudrajCopula object
##' @author Jun Yan and Ivan Kojadinovic
##'
khoudrajiCopula <- function(copula1 = indepCopula(), copula2 = indepCopula(dim=d),
                            shapes = rep(NA_real_, dim(copula1))) {

    d <- dim(copula1)

    ## if d==2, create a khoudrajiBivCopula object
    ## if d > 2
    ##          if copula1 and copula2 are explicit copulas
    ##               create a khoudrajiExplicitCopula object
    ##               (for which pdrCopula will work)
    ##          else create a khoudrajiCopula object
    ##               (for which only prCopula will work)


    ## check if copula1 and copula2 have 'exprdist' slots
    ## areNotBothExplicit <- if(is.na(match("exprdist", slotNames(copula1))) ||
    ##                       is.na(match("exprdist", slotNames(copula2))) ||
    ##                       !is.language(F1 <- copula1@exprdist$cdf) ||
    ##                       !is.language(F2 <- copula2@exprdist$cdf)) TRUE else FALSE

    ## for the moment, check if copula1 and copula2 are archmCopula, which,
    ## for the moment, implies that copula1 and copula2 have 'exprdist' slots
    areBothExplicit <- is(copula1, "archmCopula") && is(copula2, "archmCopula")

    ## d == 2 or non-explicit Khourdraji copulas
    if (d == 2 || !areBothExplicit)
        new(if (d == 2) "khoudrajiBivCopula" else "khoudrajiCopula",
            copula1 = copula1,
            copula2 = copula2,
            shapes = shapes)
    else {

        ## Explicit Khourdraji copulas

        F1 <- copula1@exprdist$cdf
        F2 <- copula2@exprdist$cdf

        ## FIXME: not characters and parse(text=), rather expressions, substitute() ...

        ## The following block handles the cdf
        ##'' @title Replace ui with (ui^shp) for each ui in a cdf expression
        ##'' @param cdf cdf expression of a copula
        ##'' @param om logical, TRUE = use 1 - shp in place of shp
        ##'' @return a string of cdf expression with (ui^shp) or (ui^(1 - shp)) in place of ui
        getcdfchar <- function(cdf, om=FALSE) {
            ## FIXME: this only works up to dim 9; e.g., u10 could be replaced with u1^shp1
            ## -----  TODO: be smarter in gsub() --- rather do *NOT* use characters at all !!
            if (d >= 10) stop("The maximum implemented dim is 9.")
            cdf <- deparse(cdf)
            for (i in 1:d) {
                ui <- paste0("u", i)
                shpi <- paste0("shp", i)
                if (om) shpi <- paste("(1 - ", shpi, ")")
                replacement <- paste("(", ui, "^", shpi, ")")
                cdf <- gsub(ui, replacement, cdf)
            }
            cdf
        }

        ## FIXME: work with expressions F1 / F2, not chars...
        cdf1 <- gsub("alpha", "c1alpha", getcdfchar(F1, om=TRUE))
        cdf2 <- gsub("alpha", "c2alpha", getcdfchar(F2, om=FALSE))
        cdf <- parse(text = c("(", cdf1, ") * (", cdf2, ")"))
        ## cdf <- substitute((F1) * (F2), list(F1 = cdf1, F2 = cdf2))
        cdf.algr <- deriv(cdf, "nothing")

        ## The following block handles pdf
        ##'' @title Get pdf expression by differentiating the cdf with D iteratively
        ##'' @param cdf Expression of cdf
        ##'' @param n dimension
        ##'' @return Expression of pdf
        pdfExpr <- function(cdf, n) {
            ## This function returns the pdf expression by differentiating the cdf
            for (i in 1:n)
                cdf <- D(cdf, paste0("u", i))
            cdf
        }

        pdf <-
            if (d <= 6)
                pdfExpr(cdf, d)
            else {
                warning("The pdf is only available for dim 6 or lower.")
                NULL
            }
        pdf.algr <- deriv(pdf, "nothing")
        exprdist <- c(cdf, pdf)
        attr(exprdist, "cdfalgr") <- cdf.algr
        attr(exprdist, "pdfalgr") <- pdf.algr

        ## The following block handles the partial derivatives of a component copula cdf
        ## needed in the density
        ## derExprs: get the derivatives of cdf of order 1 to n
        ## That is: dC/du1, d2C/(du1 du2), d3C / (du1 du2 du3), ...
        ## WARNING: This is assuming exchangeable copula so that
        ## the ordering of arguments does not matter
        derExprs <- function(cdf, n) {
            val <- rep(as.expression(cdf), n + 1) ## the first one is cdf itself
            for (i in 1:n) {
                val[i + 1] <- as.expression(D(val[i], paste0("u", i)))
                ## FIXME: why is this as.expression necessary? Martin may know.
            }
            val
        }
        derExprs1 <- derExprs(F1, d)
        derExprs2 <- derExprs(F2, d)

        new("khoudrajiExplicitCopula",
            copula1 = copula1,
            copula2 = copula2,
            shapes = shapes,
            exprdist = exprdist,
            derExprs1 = derExprs1,
            derExprs2 = derExprs2)
    }
}

## In CRAN's copula up to 0.999-14 i.e  mid-2016: --> deprecated now
asymCopula <-
asymExplicitCopula <- function(shapes, copula1, copula2) {
    .Deprecated("khoudrajiCopula(c.1, c.2, shapes) -- *reordered* arguments")
    khoudrajiCopula(copula1, copula2, shapes)
}

##################################################################################
### Methods for all Khoudraji copulas
##################################################################################

## pCopula: for all Khoudraji copulas
pKhoudrajiCopula <- function(u, copula, ...) {
    d <- dim(copula)
    tu <- if(is.matrix(u)) t(u) else matrix(u, nrow = d)
    p1 <- pCopula(t(tu^(1 - copula@shapes)), copula@copula1)
    p2 <- pCopula(t(tu^copula@shapes), copula@copula2)
    p1 * p2
}

setMethod("pCopula", signature("matrix", "khoudrajiCopula"), pKhoudrajiCopula)

## rCopula: for all Khoudraji copulas
setMethod("rCopula", signature("numeric", "khoudrajiCopula"),
          function(n, copula) {
    u <- rCopula(n, copula@copula1)
    v <- rCopula(n, copula@copula2)
    d <- dim(copula)
    x <- matrix(NA, n, d)
    ig <- KhoudFn$ig
    for (i in seq_len(d)) {
        x[,i] <- pmax(ig(u[,i], 1 - copula@shapes[i]),
                      ig(v[,i],     copula@shapes[i]))
    }
    x
})

##################################################################################
### Methods for bivariate Khoudraji copulas
##################################################################################

## dCopula: Restricted to *bivariate*  copulas
dKhoudrajiBivCopula <- function(u, copula, log = FALSE, ...) {
    a1 <- copula@shapes[1]
    a2 <- copula@shapes[2]

    ## the density can be computed only if dCdu is implemented for argument copulas
    if (!hasMethod(dCdu, class(copula@copula1)) ||
        !hasMethod(dCdu, class(copula@copula2)))
        stop("The argument copulas must both have the 'dCdu()' method implemented")

    g <- KhoudFn$g ; dgdu <- KhoudFn$dgdu
    gu1 <- cbind(g(u[,1], 1 - a1), g(u[,2], 1 - a2))
    gu2 <- cbind(g(u[,1],     a1), g(u[,2],     a2))
    dC1du <- dCdu(copula@copula1, gu1)
    dC2du <- dCdu(copula@copula2, gu2)
    ddu1_ <- dgdu(u[,1], 1 - a1)
    ddu1  <- dgdu(u[,1],     a1)
    ddu2_ <- dgdu(u[,2], 1 - a2)
    ddu2  <- dgdu(u[,2],     a2)
    part1 <- dCopula(gu1, copula@copula1) * ddu1_ * ddu2_ *
        pCopula(gu2, copula@copula2)
    part2 <- dC1du[,1] * ddu1_ * ddu2 * dC2du[,2]
    part3 <- dC1du[,2] * ddu2_ * ddu1 * dC2du[,1]
    part4 <- pCopula(gu1, copula@copula1) * dCopula(gu2, copula@copula2) * ddu2 * ddu1

    ## FIXME: use lsum() and similar to get much better numerical accuracy for log - case
    if(log)
        log(part1 + part2 + part3 + part4)
    else    part1 + part2 + part3 + part4
}

setMethod("dCopula", signature("matrix", "khoudrajiBivCopula"),
          dKhoudrajiBivCopula)

## A: Restricted to *bivariate* Khoudraji copulas
## A: Pickands dependence function only if copula1 and copula2 are extreme-value
setMethod("A", signature("khoudrajiBivCopula"), function(copula, w) {

    ## the A function can be computed only if the argument copulas are extreme-value copulas
    stopifnot(is(copula@copula1, "evCopula"), is(copula@copula2, "evCopula"))

    a1 <- copula@shapes[1];  a2 <- copula@shapes[2]
    den1 <- (1 - a1) * (1 - w) + (1 - a2) * w
    den2 <- a1 * (1 - w) + a2 * w
    t1 <- (1 - a2) * w / den1; t1[is.na(t1)] <- 1
    t2 <-      a2  * w / den2; t2[is.na(t2)] <- 1
    den1 * A(copula@copula1, t1) + den2 * A(copula@copula2, t2)
})

## dCdu: Restricted to *bivariate* Khoudraji copulas
setMethod("dCdu", signature("khoudrajiBivCopula"),
          function(copula, u, ...) {
    a1 <- copula@shapes[1]
    a2 <- copula@shapes[2]

    ## dCdu can be computed only if dCdu is implemented for argument copulas
    if (!hasMethod(dCdu, class(copula@copula1)) ||
        !hasMethod(dCdu, class(copula@copula2)))
        stop("The argument copulas must both have the 'dCdu()' method implemented")

    g <- KhoudFn$g ; dgdu <- KhoudFn$dgdu
    gu1 <- cbind(g(u[,1], 1 - a1), g(u[,2], 1 - a2))
    gu2 <- cbind(g(u[,1],     a1), g(u[,2],     a2))
    dC1du <- dCdu(copula@copula1, gu1)
    dC2du <- dCdu(copula@copula2, gu2)
    pC1gu1 <- pCopula(gu1, copula@copula1)
    pC2gu2 <- pCopula(gu2, copula@copula2)
    cbind(dgdu(u[,1], 1 - a1) * dC1du[,1] * pC2gu2 + pC1gu1 * dgdu(u[,1], a1) * dC2du[,1],
          dgdu(u[,2], 1 - a2) * dC1du[,2] * pC2gu2 + pC1gu1 * dgdu(u[,2], a2) * dC2du[,2])
})

## dCdtheta: Restricted to *bivariate* Khoudraji copulas
setMethod("dCdtheta", signature("khoudrajiBivCopula"),
          function(copula, u, ...) {
    a1 <- copula@shapes[1]
    a2 <- copula@shapes[2]
    ## dCdu can be computed only if dCdu is implemented for argument copulas
    if (!hasMethod(dCdu, class(copula@copula1)) ||
        !hasMethod(dCdu, class(copula@copula2)))
        stop("The argument copulas must both have the 'dCdu()' method implemented")
    free <- isFreeP(copula@shapes)
    g <- KhoudFn$g
    gu1 <- cbind(g(u[,1], 1 - a1), g(u[,2], 1 - a2))
    gu2 <- cbind(g(u[,1], a1), g(u[,2], a2))
    dC1du <- dCdu(copula@copula1, gu1)
    dC2du <- dCdu(copula@copula2, gu2)
    pC1gu1 <- pCopula(gu1, copula@copula1)
    pC2gu2 <- pCopula(gu2, copula@copula2)
    cbind(if(nParam(copula@copula1, freeOnly=TRUE) > 0) dCdtheta(copula@copula1, gu1) * pC2gu2,# else NULL
          if(nParam(copula@copula2, freeOnly=TRUE) > 0) pC1gu1 * dCdtheta(copula@copula2, gu2),#  "    "
          if (free[1]) -log(u[,1]) * g(u[,1], 1 - a1) * dC1du[,1] * pC2gu2 +
          pC1gu1 * log(u[,1]) * g(u[,1], a1) * dC2du[,1],
          if (free[2]) -log(u[,2]) * g(u[,2], 1 - a2) * dC1du[,2] * pC2gu2 +
          pC1gu1 * log(u[,2]) * g(u[,2], a2) * dC2du[,2])
})

##################################################################################
### dCopula method for Explicit Khoudraji copulas
##################################################################################

## This function uses the algorithmic expressions stored in the class object
dKhoudrajiExplicitCopula.algr <- function(u, copula, log=FALSE, ...) {
    dim <- dim(copula)
    stopifnot(!is.null(d <- ncol(u)), dim == d)

    for (i in 1:dim) {
        assign(paste0("u", i), u[,i])
        assign(paste0("shp", i), copula@shapes[i])
    }
    ## WARNING: assuming one parameter (or no-parameter) copula components
    c1alpha <- copula@copula1@parameters
    c2alpha <- copula@copula2@parameters
    dens <- c(eval(attr(copula@exprdist, "pdfalgr")))
    if(log) log(dens) else dens
}

setMethod("dCopula", signature("matrix",  "khoudrajiExplicitCopula"), dKhoudrajiExplicitCopula.algr)
