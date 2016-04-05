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
### Asymmetric bivariate copulas of the form C(u^(1-a)) * D(u^a)
### as in Liebscher (2008, JMVA) of which Khoudraji's device is a special case
##################################################################################

## Asymmetric copula class constructed from two d-dimensional copulas
setClass("asymCopula", contains = c("copula", "VIRTUAL"),
         representation = representation(
             copula1 = "copula",
             copula2 = "copula"
         ),
	 validity = function(object) {
    if(object@copula1@dimension != object@copula2@dimension)
        "The argument copulas are not of the same dimension"
	     else TRUE
})

##################################################################################
### Asymmetric bivariate copulas of the form
### C(u_1^{1-a_1}, u_2^{1-a_2}) * D(u_1^a_1, u_2^a_1)
##################################################################################

## Bivariate asymmetric copula class
setClass("asymBivCopula", contains = "asymCopula", ## -> calls *ITS* validity (eq. dim)
	 validity = function(object) {
    if(object@copula1@dimension != 2)
        "The asymBivCopula constituents must have dimension two"
	     else TRUE
})

## C(u_1^{1-a_1}, u_2^{1-a_2}) * D(u_1^a_1, u_2^a_1) = C(g(u, 1-a)) * D(g(u, a)
g <- function(u, a) u^a
ig <- function(u, a) u^(1 / a) # inverse
dgdu <- function(u, a) a * u^(a - 1) # derivative wrt u

##' Creates an asymBivCopula object
##'
##' @title Asymmetric bivariate copula constructor
##' @param copula1
##' @param copula2 each a bivariate exchangeable \code{\linkS4class{copula}}
##' @param shapes a numeric with two elements in [0,1]
##' @return a new "asymBivCopula" object; see above
##' @author Jun Yan and Ivan Kojadinovic
asymBivCopula <- function(copula1 = indepCopula(),
                          copula2 = indepCopula(),
                          shapes = c(1,1)) {
  new("asymBivCopula",
      dimension = copula1@dimension,
      parameters = c(copula1@parameters, copula2@parameters, shapes),
      param.names = c(if (length(copula1@parameters) > 0)
                      paste0("C1.",copula1@param.names) else character(0),
                      if (length(copula2@parameters) > 0)
                      paste0("C2.",copula2@param.names) else character(0),
                      "shape1", "shape2"),
      param.lowbnd = c(copula1@param.lowbnd, copula2@param.lowbnd, 0, 0),
      param.upbnd = c(copula1@param.upbnd, copula2@param.upbnd, 1, 1),
      copula1 = copula1,
      copula2 = copula2,
      fullname = paste("Asymmetric bivariate copula constructed from: [",
                       copula1@fullname, "] and: [", copula2@fullname, "]"))
}

## Returns shapes, copula1 and copula2 from an asymCopula object
## Not restricted to bivariate asymmetric copulas
getAsymCopulaComps <- function(object) {
    p1 <- length(object@copula1@parameters)
    p2 <- length(object@copula2@parameters)
    d <- object@dimension
    shapes <- object@parameters[(p1 + p2) + 1:d]
    copula1 <- object@copula1
    copula2 <- object@copula2
    if (p1 > 0) slot(copula1, "parameters") <- object@parameters[1:p1]
    if (p2 > 0) slot(copula2, "parameters") <- object@parameters[p1 + 1:p2]
    list(shapes = shapes, copula1 = copula1, copula2 = copula2)
}

## pCopula: Not restricted to *bivariate* asymmetric copulas
pAsymCopula <- function(u, copula) {
    tu <- t(rbind(u, deparse.level=0L))
    comps <- getAsymCopulaComps(copula)
    p1 <- pCopula(t(tu^(1 - comps$shapes)), comps$copula1)
    p2 <- pCopula(t(tu^comps$shapes), comps$copula2)
    p1 * p2
}

## rCopula: Not restricted to *bivariate* asymmetric copulas
rAsymCopula <- function(n, copula) {
    comps <- getAsymCopulaComps(copula)
    copula1 <- comps$copula1
    copula2 <- comps$copula2
    shapes <- comps$shapes
    ## Theorem 2.1, Lemma 2.1, Liebscher (2008, JMVA)
    u <- rCopula(n, copula1)
    v <- rCopula(n, copula2)
    d <- copula@dimension
    x <- matrix(NA, n, d)
    for (i in seq_len(d)) {
        x[,i] <- pmax(ig(u[,i], 1 - shapes[i]), ig(v[,i], shapes[i]))
    }
    x
}

## dCopula: Restricted to *bivariate* asymmetric copulas
dAsymBivCopula <- function(u, copula, log=FALSE, ...) {
    comps <- getAsymCopulaComps(copula)
    a1 <- comps$shapes[1]
    a2 <- comps$shapes[2]
    copula1 <- comps$copula1
    copula2 <- comps$copula2
    gu1 <- cbind(g(u[,1], 1 - a1), g(u[,2], 1 - a2))
    gu2 <- cbind(g(u[,1], a1), g(u[,2], a2))
    dC1du <- dCdu(copula1, gu1)
    dC2du <- dCdu(copula2, gu2)
    part1 <- dCopula(gu1, copula1) *
        dgdu(u[,1], 1 - a1) * dgdu(u[,2], 1 - a2) *
        pCopula(gu2, copula2)
    part2 <- dC1du[,1] * dgdu(u[,1], 1 - a1) * dgdu(u[,2], a2) * dC2du[,2]
    part3 <- dC1du[,2] * dgdu(u[,2], 1 - a2) * dgdu(u[,1], a1) * dC2du[,1]
    part4 <- pCopula(gu1, copula1) * dCopula(gu2, copula2) *
        dgdu(u[,2], a2) * dgdu(u[,1], a1)
    ## FIXME: use lsum() and similar to get much better numerical accuracy for log - case
    if(log)
        log(part1 + part2 + part3 + part4)
    else    part1 + part2 + part3 + part4
}

## A: Pickands dependence function if copula1 and copula2 are extreme-value
## Restricted to *bivariate* asymmetric copulas
AAsymCopula <- function(copula, w) {
    comps <- getAsymCopulaComps(copula)
    copula1 <- comps$copula1
    copula2 <- comps$copula2
    ## assuming copula1 and copula2 are both evCopula
    stopifnot(is(copula1, "evCopula"), is(copula2, "evCopula"))
    a1 <- comps$shapes[1];  a2 <- comps$shapes[2]
    den1 <- (1 - a1) * (1 - w) + (1 - a2) * w
    den2 <- a1 * (1 - w) + a2 * w
    t1 <- (1 - a2) * w / den1; t1 <- ifelse(is.na(t1), 1, t1)
    t2 <- a2 * w / den2; t2 <- ifelse(is.na(t2), 1, t2)
    den1 * A(copula1, t1) + den2 * A(copula2, t2)
}

## dCdu: Restricted to *bivariate* asymmetric copulas
dCduAsymBivCopula <- function(copula, u, ...) {
    comps <- getAsymCopulaComps(copula)
    a1 <- comps$shapes[1]
    a2 <- comps$shapes[2]
    copula1 <- comps$copula1
    copula2 <- comps$copula2
    gu1 <- cbind(g(u[,1], 1 - a1), g(u[,2], 1 - a2))
    gu2 <- cbind(g(u[,1], a1), g(u[,2], a2))
    dC1du <- dCdu(copula1, gu1)
    dC2du <- dCdu(copula2, gu2)
    pC1gu1 <- pCopula(gu1, copula1)
    pC2gu2 <- pCopula(gu2, copula2)
    cbind(dgdu(u[,1], 1 - a1) * dC1du[,1] * pC2gu2 +
          pC1gu1 * dgdu(u[,1], a1) * dC2du[,1],
          dgdu(u[,2], 1 - a2) * dC1du[,2] * pC2gu2 +
          pC1gu1 * dgdu(u[,2], a2) * dC2du[,2])
}

## dCdtheta: Restricted to *bivariate* asymmetric copulas
dCdthetaAsymBivCopula <- function(copula, u, ...) {
    comps <- getAsymCopulaComps(copula)
    a1 <- comps$shapes[1]
    a2 <- comps$shapes[2]
    copula1 <- comps$copula1
    copula2 <- comps$copula2
    gu1 <- cbind(g(u[,1], 1 - a1), g(u[,2], 1 - a2))
    gu2 <- cbind(g(u[,1], a1), g(u[,2], a2))
    dC1du <- dCdu(copula1, gu1)
    dC2du <- dCdu(copula2, gu2)
    pC1gu1 <- pCopula(gu1, copula1)
    pC2gu2 <- pCopula(gu2, copula2)

    if (length(copula1@parameters) > 0 && length(copula1@parameters) > 0)
        cbind(dCdtheta(copula1, gu1) * pC2gu2,
              pC1gu1 * dCdtheta(copula2, gu2),
              -log(u[,1]) * g(u[,1], 1 - a1) * dC1du[,1] * pC2gu2 +
               pC1gu1 * log(u[,1]) * g(u[,1], a1) * dC2du[,1],
              -log(u[,2]) * g(u[,2], 1 - a2) * dC1du[,2] * pC2gu2 +
               pC1gu1 * log(u[,2]) * g(u[,2], a2) * dC2du[,2])
    else if (length(copula1@parameters) == 0)
        cbind(pC1gu1 * dCdtheta(copula2, gu2),
              -log(u[,1]) * g(u[,1], 1 - a1) * dC1du[,1] * pC2gu2 +
               pC1gu1 * log(u[,1]) * g(u[,1], a1) * dC2du[,1],
              -log(u[,2]) * g(u[,2], 1 - a2) * dC1du[,2] * pC2gu2 +
               pC1gu1 * log(u[,2]) * g(u[,2], a2) * dC2du[,2])
    else cbind(-log(u[,1]) * g(u[,1], 1 - a1) * dC1du[,1] * pC2gu2 +
               pC1gu1 * log(u[,1]) * g(u[,1], a1) * dC2du[,1],
              -log(u[,2]) * g(u[,2], 1 - a2) * dC1du[,2] * pC2gu2 +
               pC1gu1 * log(u[,2]) * g(u[,2], a2) * dC2du[,2])
}

setMethod("A", signature("asymBivCopula"), AAsymCopula)

setMethod("rCopula", signature("numeric", "asymCopula"), rAsymCopula)

setMethod("pCopula", signature("numeric", "asymCopula"),pAsymCopula)
setMethod("pCopula", signature("matrix", "asymCopula"), pAsymCopula)

setMethod("dCopula", signature("numeric", "asymBivCopula"), dAsymBivCopula)
setMethod("dCopula", signature("matrix", "asymBivCopula"), dAsymBivCopula)

setMethod("dCdu", signature("asymBivCopula"), dCduAsymBivCopula)
setMethod("dCdtheta", signature("asymBivCopula"), dCdthetaAsymBivCopula)


################################################################################
### Asymmetric *explicit* copulas
### Explicit refers to the fact that the two d-dimensional copulas
### have explicit cdfs and pdfs
################################################################################

setClass("asymExplicitCopula", contains = "asymCopula",
         representation = representation(
           exprdist = "expression",
           derExprs1 = "expression",
           derExprs2 = "expression"
           )
         ## validity = function(object) {
         ##     ## TODO: check exprdist, derExprs[12]
         ## },
         )

##' Creates an asymExplicitCopula object
##'
##' @title Asymmetric explicit copula constructor
##' @param copula1 a d-dimensional copula
##' @param copula2 a d-dimensional copula
##' @param shapes a numeric with elements in [0,1]
##' @return a new "asymExplicitCopula" object; see above
##' @author Jun Yan and Ivan Kojadinovic
asymExplicitCopula <- function(copula1, copula2, shapes) {
    stopifnot(copula2@dimension == (d <- copula1@dimension),
              length(shapes) == d)
    ## check if have 'exprdist' slots [otherwise "bad" error msg below]:
    if(is.na(match("exprdist", slotNames(copula1))))
        stop(gettextf("Not yet implemented for '%s' of class \"%s\"",
                      "copula1", class(copula1)), domain = NA)
    if(is.na(match("exprdist", slotNames(copula2))))
        stop(gettextf("Not yet implemented for '%s' of class \"%s\"",
                      "copula2", class(copula2)), domain = NA)
    stopifnot(is.language(F1 <- copula1@exprdist$cdf),
              is.language(F2 <- copula2@exprdist$cdf))

    ## FIXME: not characters and parse(text=), rather expressions, substitute() ...

    ## cdf
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
    cdf1 <- getcdfchar(F1, om=TRUE)
    cdf2 <- getcdfchar(F2, om=FALSE)
    cdf <- parse(text = c("(", cdf1, ") * (", cdf2, ")"))
    ## cdf <- substitute((F1) * (F2),
    ##                   list(F1 = cdf1, F2 = cdf2))

    ## pdf
    pdfExpr <- function(cdf, n) {
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
    derExprs <- function(cdf, n) {
        val <- as.expression(cdf)
        for (i in 1:n) {
            val <- c(val, D(val[i], paste0("u", i)))
        }
        val
    }
    derExprs1 <- derExprs(F1, d)
    derExprs2 <- derExprs(F2, d)
    shapes.names <- paste0("shape", 1:d)
    new("asymExplicitCopula",
        dimension = d,
        parameters = c(copula1@parameters, copula2@parameters, shapes),
        param.names = c(copula1@param.names, copula2@param.names, shapes.names),
        param.lowbnd = c(copula1@param.lowbnd, copula2@param.lowbnd, rep(0, d)),
        param.upbnd  = c(copula1@param.upbnd,  copula2@param.upbnd,  rep(1, d)),
        copula1 = copula1,
        copula2 = copula2,
        exprdist = c(cdf=cdf, pdf=pdf),
        derExprs1 = derExprs1, derExprs2 = derExprs2,
        fullname = "Asymmetric Explicit Copula")
}

getPowerSet <- function(d) {
  TF <- matrix(c(TRUE, FALSE), 2, d)
  as.matrix(expand.grid(as.list(as.data.frame(TF))))
}

densDers <- function(idx, u, dg, copula, derExprs) {
    ## assuming exchangeable copula1 and copula2
    dorder <- sum(idx)
    alpha <- copula@parameters[1] # possibly needed in 'derExprs' below
    d <- copula@dimension
    newidx <- c((1:d)[idx], (1:d)[!idx])
    u <- u[, newidx]
    for (i in 1:d) assign(paste0("u", i), u[,i])
    dgu <- if(sum(idx) == 0) 1 else apply(dg[,idx,drop=FALSE], 1, prod)
    c(eval(derExprs[dorder + 1])) * dgu
}

dAsymExplicitCopula <- function(u, copula, log=FALSE, ...) {
    u <- as.matrix(u)
    comps <- getAsymCopulaComps(copula)
    a <- comps$shapes
    tu <- t(u)
    u1 <- t(tu^(1 - a))
    u2 <- t(tu^a)
    dg1 <- t((1 - a) * tu^(-a))
    dg2 <- t(a * tu^(a - 1))
    d <- copula@dimension
    powerSet <- getPowerSet(d)
    dens <- 0
    for (i in 1:nrow(powerSet)) {
        idx1 <- c(powerSet[i,])
        idx2 <- c(!powerSet[i,])
        part1 <- densDers(idx1, u1, dg1, copula@copula1, copula@derExprs1)
        part2 <- densDers(idx2, u2, dg2, copula@copula2, copula@derExprs2)
        dens <- dens + part1 * part2
        ## print(part1); print(part2)
    }
    if(log) log(dens) else dens
}

setMethod("dCopula", signature("numeric", "asymExplicitCopula"), dAsymExplicitCopula)
setMethod("dCopula", signature("matrix", "asymExplicitCopula"), dAsymExplicitCopula)
