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
### Partial derivatives of the CDF wrt arguments
##################################################################################

setGeneric("dCdu", function(copula, u, ...) standardGeneric("dCdu"))

## Function returning a list of control arguments for grad() / jacobian()
## in package numDeriv
gradControl <- function(eps = 1e-4, d = 1e-4,
                        zero.tol = sqrt(.Machine$double.eps / 7e-7),
                        r = 6, v = 2, show.details = FALSE) {
    list(eps = eps, d = d, zero.tol = zero.tol, r = r, v = v,
         show.details = show.details)
}

## Returns 'side' vector for grad() / jacobian() when computing
## partial derivatives wrt x in (lb, ub) or [lb, ud)
## d comes from gradControl
sides <- function(x, dimx, d, lb, ub) {
    res <- rep(NA, dimx) # two-sided derivatives by default
    res[x - lb < d] <- 1 # right derivative
    res[ub - x < d] <- -1 # left derivative
    res
}

## Basic implementation based on numerical differentiation
dCduCopulaNum <- function(copula, u, method.args = gradControl(d = 1e-1), ...) {
    warning("Function dCdu() not implemented for this copula: numerical differentiation used")
    logC <- function(x) log(pCopula(x, copula))
    dim <- copula@dimension
    res <- t(apply(u, 1, function(u.)
        grad(logC, u., side = sides(u., dim, method.args$d, 0, 1),
             method.args = method.args))) * pCopula(u, copula)
    res[res < 0] <- 0 # because 0 < dCdu < 1 uniformly
    res[res > 1] <- 1
    res
}

dCduArchmCopula <- function(copula, u, ...) {
    ## TODO: For ACs, the following is better than 'fixed formulas' => use it
    ## require(copula)
    ## n <- 10
    ## d <- 4
    ## family <- "Gumbel"
    ## th <- 2
    ## u <- matrix(runif(n*d), ncol=d)
    ## cop <- onacopulaL(family, nacList=list(th, seq_len(d)))
    ## iPsi.u <- cop@copula@iPsi(u, theta=th)
    iPsi.u <- iPsi(copula, u)
    d <- copula@dimension
    j <- ceiling(d/2)
    ## NEED function absdPsi to make below work for any Archimedean copula in the same way that iPsi exists
    ## dCdu <- sapply(seq_len(d), function(j) exp(cop@copula@absdPsi(rowSums(iPsi.u), theta=th, degree=1, log=TRUE) - cop@copula@absdPsi(iPsi.u[,j], theta=th, degree=1, log=TRUE)))
}

## Warning: This function assumes symmetry in u
dCduExplicitCopula <- function(copula, u, ...)
{
    p <- copula@dimension
    mat <- matrix(NA_real_, nrow(u),p)
    algNm <- paste(class(copula)[1], "cdfDerWrtArg.algr", sep=".")
    if(exists(algNm)) {
	alpha <- copula@parameters # typically used in 'eval(der.cdf.u, *)' below
        der.cdf.u <- get(algNm)[p]
        unames0 <- paste0("u",1:p)
        for (j in 1:p)
        {
            unames <- unames0; unames[1] <- unames0[j]; unames[j] <- unames0[1]
            colnames(u) <- unames
            mat[,j] <- eval(der.cdf.u, data.frame(u))
        }
    } else warning("there is no formula for dCdu() for this copula")
    mat
}


## this function is used for Khoudraji's device
dCduIndepCopula <- function(copula, u, ...) {
  p <- copula@dimension
  mat <- matrix(0, nrow(u), p)
  for (j in 1:p) {
    mat[,j] <- apply(u[, -j, drop=FALSE], 1, prod)
  }
  mat
}

setMethod("dCdu", signature("Copula"), dCduCopulaNum)
setMethod("dCdu", signature("archmCopula"), dCduExplicitCopula)
setMethod("dCdu", signature("plackettCopula"), dCduExplicitCopula)
setMethod("dCdu", signature("evCopula"), dCduExplicitCopula)
setMethod("dCdu", signature("gumbelCopula"), dCduExplicitCopula)
setMethod("dCdu", signature("indepCopula"), dCduIndepCopula)

dCduEllipCopula <- function(copula, u, ...)
{
    p <- copula@dimension
    sigma <- getSigma(copula)

    ## quantile transformation
    v <- switch(class(copula),
		"normalCopula" = qnorm(u),
		"tCopula" = {
		    df <- copula@df
		    qt(u, df=df)
		},
		stop("not implemented for class ", class(copula)))
    n <- nrow(u)
    mat <- matrix(0,n,p)

    for (j in 1:p)
    {
	s <- sigma[-j,-j] - sigma[-j,j,drop=FALSE] %*% sigma[j,-j,drop=FALSE]

	switch(class(copula),
               "normalCopula" =
           {
               if (p == 2) {
                   rho <- copula@parameters
                   mat[,j] <- pnorm(v[,-j], rho * v[,j], sqrt(1 - rho^2))
               }
               else
                   for (i in 1:n)
                       mat[i,j] <- pmvnorm(lower = rep(-Inf, p - 1), upper = v[i,-j],
                                           mean = v[i,j] * sigma[-j,j],
                                           sigma = drop(s))
           },
               "tCopula" =
           {
               if (p == 2) {
                   rho <- copula@parameters
                   mat[,j] <-  pt(sqrt((df+1)/(df+v[,j]^2)) / sqrt(1 - rho^2)
                                  * (v[,-j] - rho * v[,j]), df=df+1)
               }
               else {
                   if(df != as.integer(df))
                       stop("'df' is not integer; therefore, dCdu() cannot be computed yet")
                   for (i in 1:n)
                       mat[i,j] <- pmvt(lower = rep(-Inf, p - 1),
                                        upper = drop(sqrt((df+1)/(df+v[i,j]^2)) *
                                        (v[i,-j] - v[i,j] * sigma[-j,j])),
                                        sigma = s, df = df + 1)
               }

           })
    }
    mat
}

setMethod("dCdu", signature("ellipCopula"), dCduEllipCopula)

##################################################################################
### Plackett formula for elliptical copulas
##################################################################################

setGeneric("plackettFormulaDim2", function(copula, x) standardGeneric("plackettFormulaDim2"))

plackettFormulaDim2NormalCopula <- function(copula, x)
  {
    rho <- copula@parameters
    ir2 <- 1 - rho^2
    as.matrix(exp(-(x[,1]^2 + x[,2]^2 - 2 * rho * x[,1] * x[,2]) /
                  (2 * ir2)) /
              (2 * pi * sqrt(ir2)))
  }

setMethod("plackettFormulaDim2", signature("normalCopula"), plackettFormulaDim2NormalCopula)

plackettFormulaDim2TCopula <- function(copula, x)
  {
    rho <- copula@parameters
    ir2 <- 1 - rho^2
    df <- copula@df
    as.matrix((1 + (x[,1]^2 + x[,2]^2 - 2 * rho * x[,1] * x[,2]) /
               (df * ir2))^(-df / 2) / (2 * pi * sqrt(ir2)))
  }

setMethod("plackettFormulaDim2", signature("tCopula"), plackettFormulaDim2TCopula)

setGeneric("plackettFormula",  function(copula, p, rho, s, m, x, i, j) standardGeneric("plackettFormula"))

plackettFormulaNormalCopula <- function(copula, p, rho, s, m, x, i, j)
{
    exp(-(x[i]^2 + x[j]^2 - 2 * rho * x[i] * x[j]) /
               (2 * (1 - rho^2))) / (2 * pi * sqrt(1 - rho^2)) *
           (if (p == 3) pnorm(drop((x[-c(i,j)] - m %*% x[c(i,j)])/sqrt(s)))
           else pmvnorm(lower = rep(-Inf, p - 2),
                        upper = drop(x[-c(i,j)] - m %*% x[c(i,j)]),
                        sigma = s))
}

setMethod("plackettFormula", signature("normalCopula"), plackettFormulaNormalCopula)

plackettFormulaTCopula <- function(copula, p, rho, s, m, x, i, j)
{
    stopifnot(p >= 3)
    df <- copula@df
    if(df != as.integer(df) && p > 3)
	stop("'df' is not integer; therefore, plackettFormula() cannot be computed yet")
    term <- 1 + (x[i]^2 + x[j]^2 - 2 * rho * x[i] * x[j]) / (df * (1 - rho^2))
    ## return:
    term^(-df / 2) / (2 * pi * sqrt(1 - rho^2)) *
	if (p == 3) pt(drop((x[-c(i,j)] - m %*% x[c(i,j)]) / sqrt(term * s)), df= df)
	else pmvt(df = df, lower = rep(-Inf, p - 2),
		  upper = drop((x[-c(i,j)] - m %*% x[c(i,j)]) / sqrt(term)),
		  sigma = s)
}

setMethod("plackettFormula", signature("tCopula"), plackettFormulaTCopula)

##################################################################################
### Partial derivatives of the CDF wrt parameters
##################################################################################

setGeneric("dCdtheta", function(copula, u, ...) standardGeneric("dCdtheta"))

## Basic implementation based on numerical differentiation
dCdthetaCopulaNum <- function(copula, u, method.args = gradControl(d = 1e-1), ...) {
    warning("Function dCdtheta() not implemented for this copula: numerical differentiation used")
    logC <- function(theta) {
        copula@parameters <- theta
        log(pCopula(u, copula))
    }
    theta <- copula@parameters
    p <- length(theta)
    lb <- copula@param.lowbnd
    ub <- copula@param.lowbnd
    jacobian(logC, theta, side = sides(theta, p, method.args$d, lb, ub),
             method.args = method.args) * pCopula(u, copula)
}

dCdthetaExplicitCopula <- function(copula, u, ...)
{
    p <- copula@dimension
    algNm <- paste(class(copula)[1], "cdfDerWrtPar.algr", sep=".")
    if(exists(algNm)) {
        alpha <- copula@parameters # typically used in 'eval(*)' below
        der.cdf.alpha <- get(algNm)[p]
        colnames(u) <- paste0("u", 1:p)
        as.matrix(eval(der.cdf.alpha, data.frame(u)))
    } else {
        warning("there is no formula for dCdtheta() for this copula")
        matrix(NA_real_, nrow(u),p)
    }
}

dCdthetaEvCopula <- function(copula, u, ...) {
  loguv <- log(u[,1], u[,2])
  w <- log(u[,2]) / loguv
  ## return  der.cdf.alpha
  as.matrix(pCopula(u, copula) * loguv * dAdtheta(copula, w))
}

setMethod("dCdtheta", signature("Copula"), dCdthetaCopulaNum)
setMethod("dCdtheta", signature("archmCopula"), dCdthetaExplicitCopula)
setMethod("dCdtheta", signature("plackettCopula"), dCdthetaExplicitCopula)
setMethod("dCdtheta", signature("evCopula"), dCdthetaExplicitCopula)
setMethod("dCdtheta", signature("gumbelCopula"), dCdthetaExplicitCopula)

dCdthetaEllipCopula <- function(copula, u, ...)
{
    p <- copula@dimension

    ## quantile transformation
    v <- switch(class(copula),
		"normalCopula" = qnorm(u),
		"tCopula" = qt(u, df = copula@df),
		stop("not implemented for class ", class(copula)))

    if (p == 2)
        plackettFormulaDim2(copula, v)
    else
    {
        n <- nrow(u)
        sigma <- getSigma(copula)

        if (copula@dispstr %in% c("ex","ar1")) ## exchangeable or ar1
        {
            rho <- copula@parameters
            r <- matrix(c(1,-rho,-rho,1),2,2)/(1 - rho^2)
            m <- sigma[-c(1,2),c(1,2)] %*% r
            s <- sigma[-c(1,2),-c(1,2)] - sigma[-c(1,2),c(1,2)] %*% r %*% sigma[c(1,2),-c(1,2)]

            mat <- matrix(0,n,1)

            if (copula@dispstr == "ex") ## exchangeable
              for (k in 1:n)
                for (j in 1:(p-1))
                  for (i in (j+1):p)
                    mat[k,1] <- mat[k,1] + plackettFormula(copula, p, rho, s, m, v[k,], i, j)
            else ## ar1
              for (k in 1:n)
                for (j in 1:(p-1))
                  for (i in (j+1):p)
                    mat[k,1] <- mat[k,1] + (i - j) * rho ^ (i - j - 1) *
				  plackettFormula(copula, p, rho, s, m, v[k,], i, j)

            mat
        }
	else # unstructured or toeplitz or ...
        {
            mat <- matrix(0,n,p*(p-1)/2)
            l <- 1
            for (j in 1:(p-1))
                for (i in (j+1):p)
                {
                    rho <- sigma[i,j]
                    r <- matrix(c(1,-rho,-rho,1),2,2)/(1 - rho^2)
                    m <- sigma[-c(i,j),c(i,j)] %*% r
                    s <- sigma[-c(i,j),-c(i,j)] - m %*% sigma[c(i,j),-c(i,j)]

                    for (k in 1:n)
                        mat[k,l] <- plackettFormula(copula, p, rho, s, m, v[k,], i, j)
                    l <- l + 1
                }
            if (copula@dispstr == "un") ## unstructured
                mat
            else if (copula@dispstr == "toep") {
                coef <- matrix(0,p*(p-1)/2,p-1)
                for (k in 1:(p-1))
                {
                    m <- row(sigma) == col(sigma) + k
                    coef[,k] <- P2p(m)
                }
                mat %*% coef
            }
            else stop("Not implemented yet for the dispersion structure ", copula@dispstr)
        }
    }
}

setMethod("dCdtheta", signature("ellipCopula"), dCdthetaEllipCopula)

##################################################################################
### Partial derivatives of the log PDF wrt arguments
##################################################################################

setGeneric("dlogcdu", function(copula, u, ...) standardGeneric("dlogcdu"))

## Basic implementation based on numerical differentiation
dlogcduCopulaNum <- function(copula, u, method.args = gradControl(d = 1e-5), ...) {
    warning("Function dlogcdu() not implemented for this copula: numerical differentiation used")
    logc <- function(x) dCopula(x, copula, log = TRUE)
    dim <- copula@dimension
    t(apply(u, 1, function(u.)
        grad(logc, u., side = sides(u., dim, method.args$d, 0, 1),
             method.args = method.args)))
}

dlogcduExplicitCopula <- function(copula, u, ...)
{
    p <- copula@dimension
    algNm <- paste(class(copula)[1], "pdfDerWrtArg.algr", sep=".")
    mat <- matrix(NA_real_, nrow(u), p)
    if(exists(algNm)) {
        der.pdf.u <- get(algNm)[p]
	alpha <- copula@parameters # typically used in eval()
        unames0 <- paste0("u", 1:p)
        for (j in 1:p)
        {
            unames <- unames0; unames[1] <- unames0[j]; unames[j] <- unames0[1]
            colnames(u) <- unames
            mat[,j] <- eval(der.pdf.u, data.frame(u))
        }
    } else warning("there is no formula for dlogcdu() for this copula")
    mat / dCopula(u, copula)
}

setMethod("dlogcdu", signature("Copula"), dlogcduCopulaNum)
setMethod("dlogcdu", signature("archmCopula"), dlogcduExplicitCopula)
setMethod("dlogcdu", signature("plackettCopula"), dlogcduExplicitCopula)
setMethod("dlogcdu", signature("evCopula"), dlogcduExplicitCopula)
setMethod("dlogcdu", signature("gumbelCopula"), dlogcduExplicitCopula)

dlogcduNormalCopula <- function(copula, u, ...)
{
    v <- qnorm(u)
    (- v %*% solve(getSigma(copula)) + v) / dnorm(v)
}

setMethod("dlogcdu", signature("normalCopula"), dlogcduNormalCopula)

dlogcduTCopula <- function(copula, u, ...)
{
    df <- copula@df
    v <- qt(u,df=df)
    w <- dt(v,df=df)
    m <- v %*% solve(getSigma(copula))
    - (df + copula@dimension) * m / ((df + rowSums(m * v)) * w) +
        (df + 1) * v / ((df +  v^2) * w)
}

setMethod("dlogcdu", signature("tCopula"), dlogcduTCopula)

##################################################################################
### Partial derivatives of the log PDF wrt parameters
##################################################################################

setGeneric("dlogcdtheta", function(copula, u, ...) standardGeneric("dlogcdtheta"))

## Basic implementation based on numerical differentiation
dlogcdthetaCopulaNum <- function(copula, u, method.args = gradControl(d = 1e-5), ...) {
    warning("Function dlogcdtheta() not implemented for this copula: numerical differentiation used")
    logc <- function(theta) {
        copula@parameters <- theta
        dCopula(u, copula, log = TRUE)
    }
    theta <- copula@parameters
    p <- length(theta)
    lb <- copula@param.lowbnd
    ub <- copula@param.lowbnd
    jacobian(logc, theta, side = sides(theta, p, method.args$d, lb, ub),
             method.args = method.args)
}

dlogcdthetaExplicitCopula <- function(copula, u, ...)
{
    p <- copula@dimension
    algNm <- paste(class(copula)[1], "pdfDerWrtPar.algr", sep=".")
    if(exists(algNm) && !is.null((der.pdf.alpha <- get(algNm)[p])[[1]])) {
	alpha <- copula@parameters # typically used in val(.)
        colnames(u) <- paste0("u", 1:p)
        as.matrix(eval(der.pdf.alpha, data.frame(u))) / dCopula(u, copula)
    } else {
        warning("There is no formula for dlogcdtheta() for this copula")
        matrix(NA_real_, nrow(u), p)
    }
}

setMethod("dlogcdtheta", signature("Copula"), dlogcdthetaCopulaNum)
setMethod("dlogcdtheta", signature("archmCopula"), dlogcdthetaExplicitCopula)
setMethod("dlogcdtheta", signature("plackettCopula"), dlogcdthetaExplicitCopula)
setMethod("dlogcdtheta", signature("evCopula"), dlogcdthetaExplicitCopula)

dlogcdthetaEllipCopula <- function(copula, u, ...)
{
    p <- copula@dimension

    ## quantile transformation
    v <- switch(clc <- class(copula),
		"normalCopula" = qnorm(u),
		"tCopula" = {
		    df <- copula@df
		    qt(u, df=df)
		},
		## else:
		stop("not implemented"))
    if (p == 2)
    {
	rho <- copula@parameters
	ir2 <- 1 - rho^2
        sv2 <- rowSums(v^2) # == v[,1]^2 + v[,2]^2
	as.matrix(switch(clc,
			 "normalCopula" =
			 (rho * ir2 - rho * sv2 + (rho^2 + 1) * v[,1] * v[,2])/ir2^2,
			 "tCopula" =
			 (1 + df) * rho / -ir2 + (2 + df) * (df * rho + v[,1] * v[,2])
			 / (df * ir2 + sv2 - 2 * rho * v[,1] * v[,2])))

    } else { ##  p >= 3
        n <- nrow(u)
        sigma <- getSigma(copula)
        detsig <- det(sigma)
        invsig <- solve(sigma)

        if (copula@dispstr %in% c("ex","ar1")) ## exchangeable or ar1
        {
            rho <- copula@parameters

            dersig <- matrix(1,p,p)
            if (copula@dispstr == "ex") ## ex
                diag(dersig) <- 0
            else ## ar1
                for (i in 1:p)
                    for (j in 1:p) {
                        ij <- abs(i - j)
                        dersig[i,j] <- ij * rho^(ij - 1)
                    }

            ## MM:  sum(diag(A %*% B)) == sum(A * t(B)) .. and B=dersig is symmetric here
            ## derdetsig <- detsig * sum(diag(invsig %*% dersig))
            derdetsig <- detsig *  sum(invsig * dersig)
            derinvsig <- - invsig %*% dersig %*% invsig
            firstterm <- derdetsig/detsig

	    mat <-
		switch(clc,
		       "normalCopula" =
		       - (firstterm + rowSums((v %*% derinvsig) * v))/2,
		       "tCopula" =
		       - (firstterm + (df + p) * rowSums((v %*% derinvsig) * v)
			  / (df +  rowSums((v %*% invsig) * v)) ) / 2)
            as.matrix(mat)
        }
	else # unstructured or toeplitz or ...
	  {
            mat <- matrix(0,n,p*(p-1)/2)
            l <- 1
            for (j in 1:(p-1))
                for (i in (j+1):p)
                {
                    derdetsig <- 2 * det(sigma[-i,-j,drop=FALSE]) * (-1)^(i+j)
                    derinvsig <- - invsig[,i] %*% t(invsig[,j]) - invsig[,j] %*% t(invsig[,i])
                    firstterm <- derdetsig/detsig

		    mat[,l] <-
			switch(clc,
			       "normalCopula" =
			       - (firstterm + rowSums((v %*% derinvsig) * v))/2,
			       "tCopula" =
			       - (firstterm + (df + p) * rowSums((v %*% derinvsig) * v)
				  / (df +  rowSums((v %*% invsig) * v)) ) / 2)
                    l <- l + 1
                }
            if (copula@dispstr == "un") ## unstructured
                mat
            else if (copula@dispstr == "toep") { ## toeplitz: p-1 parameters
                coef <- matrix(0,p*(p-1)/2,p-1)
                for (k in 1:(p-1))
                {
                    m <- row(sigma) == col(sigma) + k
                    coef[,k] <- P2p(m)
                }
                mat %*% coef
            }
            else stop("Not implemented yet for the dispersion structure ", copula@dispstr)
        }
    } ## p >= 3
}

setMethod("dlogcdtheta", signature("ellipCopula"), dlogcdthetaEllipCopula)
