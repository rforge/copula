##' A list of parCopula's
setClass("parClist",  contains = "list",
	 validity = function(object) {
	     if(!length(object))
		 return("empty parCopula lists are not valid")
	     cls  <- vapply(object, class, "")
	     is.pC <- vapply(cls, extends, NA, class2 = "parCopula")
	     if(!all(is.pC))
		 return(paste("components", paste(which(!is.PC), collapse=", "),
			      "do not inherit from \"parCopula\""))
	     dims <- vapply(object, dim, 1)
	     if(any(ne.d <- dims[1] != dims))
		 return(paste("copula dimensions differ", dims[1], "!=",
			      dims[match(FALSE, ne.d)]))
	     ## else
	     TRUE
	 })

##' Mixture Weights (not yet exported)
setClass("mixWeights", contains = "numeric",
	 validity = function(object) {
	     if(any(object < 0))
		 return("mixture weights must not be negative")
	     s <- sum(object)
	     if(abs(s - 1) > 64*.Machine$double.eps)
		 return(paste("mixture weights do not add to one, but to 1 -",
			      formatC(1-s)))
	     TRUE
	 })

##' A Mixture of Copulas
setClass("mixCopula", contains = "parCopula",
	 slots = c("w" = "mixWeights", "cops" = "parClist"),
	 validity = function(object) {
	     if((n1 <- length(object@w)) != (n2 <- length(object@cops)))
		 return(paste("length of weights and copula list differ:",
			      paste(n1, "!=", n2)))
	     TRUE
	 })

##' A Mixture of Explicit Copulas
setClass("mixExplicitCopula", contains = "mixCopula",
         slot = c("exprdist" = "expression")
         )


mixCopula <- function(coplist, w = NULL) {
    stopifnot(is.list(coplist))
    if((m <- length(coplist)) == 0)
	stop("a mixCopula must have at least one copula component")
    if (is.null(w)) # default: equal weights
	w <- rep.int(1/m, m)
    allExplicit <- all(vapply(coplist, isExplicit, NA))

    ## now the validity methods check:
    if (!allExplicit)
        new("mixCopula",
            w = as(w, "mixWeights"),
            cops = as(coplist, "parClist"))
    else { ## mixExplicitCopula
        mixCdf <- mixPdf <- parse(text = paste0("pcdf", 1L:m, " * ", "p", 1L:m,
                                                collapse = " + "))
        exprdist <- c(cdf = mixCdf, pdf = mixPdf)
        cdfL <- lapply(coplist, function(x) x@exprdist$cdf)
        pdfL <- lapply(coplist, function(x) x@exprdist$pdf)
        ListExpr <- function(nm1, nm2) # ((uglyness in one place))
            parse(text = paste0("list(", paste0(nm1, "= quote(",nm2,")", collapse= ", "), ")"))
        for (i in 1:m) {
            ## original 6 basic families have alpha in expressions
            cdfL[[i]] <- do.call(substitute, list(cdfL[[i]], list(alpha = quote(param))))
            pdfL[[i]] <- do.call(substitute, list(pdfL[[i]], list(alpha = quote(param))))
            ## rename the parameters with a prefix of 'm<i>.'
            oldParNames <- names(getParam(coplist[[i]], freeOnly=FALSE, named=TRUE))
            npar <- length(oldParNames)
            if (npar > 0) {
                prefix <- paste0("m", i, ".")
                newParNames <- paste0(prefix, oldParNames)
                rep.l <- ListExpr(oldParNames, newParNames)
                cdfL[[i]] <- do.call(substitute, list(cdfL[[i]], eval(rep.l)))
                pdfL[[i]] <- do.call(substitute, list(pdfL[[i]], eval(rep.l)))
            }
        }
        cdfL <- as.expression(cdfL)
        pdfL <- as.expression(pdfL)
        pcdfs <- paste0("pcdf", 1:m)
        cdf.repl <- ListExpr(pcdfs, cdfL)
        pdf.repl <- ListExpr(pcdfs, pdfL)
        ## why this does not work? what happened when they were put together with c?
        ## mixCdf <- do.call(substitute, list(mixCdf, eval(cdf.repl)))
        ## mixPdf <- do.call(substitute, list(mixPdf, eval(pdf.repl)))
        mixCdf <- do.call(substitute, list(exprdist$cdf, eval(cdf.repl)))
        mixPdf <- do.call(substitute, list(exprdist$pdf, eval(pdf.repl)))

        cdf <- as.expression(mixCdf)
        cdf.algr <- deriv(cdf, "nothing")
        pdf <- as.expression(mixPdf)
        pdf.algr <- deriv(pdf, "nothing")
        exprdist <- c(cdf = cdf, pdf = pdf)
        attr(exprdist, "cdfalgr") <- cdf.algr
        attr(exprdist, "pdfalgr") <- pdf.algr

        new("mixExplicitCopula",
            w = as(w, "mixWeights"),
            cops = as(coplist, "parClist"),
            exprdist = exprdist)
    }
}

setMethod("dim", signature("mixCopula"), function(x) dim(x@cops[[1]]))

## get parameters
setMethod("getParam", "mixCopula",
	  function(copula, freeOnly = TRUE, attr = FALSE, named = attr) {
	      d <- dim(copula)
	      parC <- lapply(copula@cops, getParam,
			     freeOnly=freeOnly, attr=TRUE, named=named)
	      m <- length(w <- copula@w)
              if(named) {
                  parC <- sapply(1:m, function(i) {
                      if (length(parC[[i]]) > 0) setNames(parC[[i]], paste0("m", i, ".", names(parC[[i]])))
                  })
                  attr(w, "names") <- paste0("p", 1L:m)
              }

	      ## FIXME re-parametrize 'w' a la nor1mix:: (??)
	      if(attr) { # more structured result
		  ## need attr(*, "param.(up|low)bnd") for for loglikCopula()
		  lowb <- lapply(parC, attr, "param.lowbnd")
		  uppb <- lapply(parC, attr, "param.upbnd")
		  structure(c(unlist(parC),                    w),
			    param.lowbnd = c(unlist(lowb), rep(0, m)),
			    param.upbnd  = c(unlist(uppb), rep(1, m)))
	      } else {
		  c(unlist(parC), w)
	      }
	  })


setMethod("freeParam<-", signature("mixCopula", "numeric"),
	  function(copula, value) {
	      cops <- copula@cops
	      m <- length(cops)
	      nj <- vapply(cops, nParam, 1, freeOnly=TRUE)
	      ## FIXME re-parametrize a la nor1mix::nM2par / .par2nM
	      ## ----- i.e. value would only contain  qlogis(w[-1])  !!
	      nw <- length(iF.w <- isFreeP(w <- copula@w))
	      if (sum(nj) + nw != length(value))
		  stop(sprintf(
		  "length(value) = %d  !=  %d, the number of free parameters",
		  length(value), sum(nj) + nw))
	      n. <- 0L
	      for(j in seq_len(m))
		  if ((n.j <- nj[j])) {
		      freeParam(cops[[j]]) <- value[n.+ seq_len(n.j)]
		      n. <- n. + n.j
		  }
	      if(n.)
		  copula@cops <- cops
	      if(nw)
		  copula@w[iF.w] <- value[n. + seq_len(nw)]
	      copula
	  })

## logical indicating which parameters are free
setMethod("isFree", signature("mixCopula"), function(copula)
    c(unlist(lapply(copula@cops, isFree)), isFreeP(copula@w))) # FIXME reparametrize 'w'

## number of (free / all) parameters :
setMethod("nParam", signature("mixCopula"), function(copula, freeOnly=FALSE)
    sum(vapply(copula@cops, nParam, 1, freeOnly=freeOnly)) +
    (if(freeOnly) nFree else length)(copula@w)) # FIXME reparametrize 'w'



setMethod(describeCop, c("mixCopula", "character"), function(x, kind) {
    m <- length(x@w)
    c1 <- paste("mixCopula from", m, "components")
    if(kind == "very short")
        return(c1)
    ## else
    dC <- vapply(x@cops, describeCop, "", kind="very short")
    paste0(c1, "\n",
	   dputNamed(if(kind == "short") unname(vapply(dC, abbreviate, "")) else dC),
	   "  with weights:\n", dputNamed(x@w))
})


##' The  C() function :
pMixCopula <- function(u, copula, ...) {
    as.vector(
	vapply(copula@cops, pCopula, FUN.VALUE=numeric(nrow(u)), u=u, ...)
	%*%
	copula@w)
}

setMethod("pCopula", signature("matrix",  "mixCopula"), pMixCopula)

##' The  c() function :
dMixCopula <- function(u, copula, log = FALSE, ...) {
    as.vector(
	vapply(copula@cops, dCopula, FUN.VALUE=numeric(nrow(u)), u=u, log=log, ...)
	%*%
	copula@w)
}

setMethod("dCopula", signature("matrix",  "mixCopula"), dMixCopula)

## Random Number Generation
setMethod("rCopula", signature("numeric", "mixCopula"),
	  function(n, copula) {
	      m <- length(w <- copula@w)
	      if(n == 1) {
		  j <- sample(m, size = 1, prob = w)
		  rCopula(1, copula@cops[[j]])
	      } else {
		  nj <- as.vector(rmultinom(n=1, size = n, prob = w))
		  ## sample nj from j-th copula
		  U <- lapply(seq(along=nj),
			      function(j) rCopula(nj[j], copula@cops[[j]]))
		  ## bind rows, and permute finally:
		  do.call(rbind, U)[sample.int(n), ]
	      }
})

setMethod("lambda", "mixCopula", function(copula, ...)
    setNames(c(vapply(copula@cops, lambda, numeric(2)) %*% copula@w),
             c("lower", "upper")))


setMethod("rho", "mixCopula", function(copula, ...)
    c(vapply(copula@cops, rho, 1.1) %*% copula@w))

