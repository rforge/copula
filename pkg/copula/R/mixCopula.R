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

##' Mixture Weights
setClass("mixWeights",  contains = "numeric",
	 validity = function(object) {
	     if(any(object < 0))
		 return("mixture weights must not be negative")
	     s <- sum(object)
	     if(abs(s - 1) > 10*.Machine$double.eps)
		 return("mixture weights add to one")
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


mixCopula <- function(coplist, w = NULL) {
    stopifnot(is.list(coplist))
    if((m <- length(coplist)) == 0)
	stop("a mixCopula must have at least one copula component")
    if (is.null(w)) # default: equal weights
	w <- rep.int(1/m, m)
    ## now the validity methods check:
    new("mixCopula",
	w = as(w, "mixWeights"),
	cops = as(coplist, "parClist"))
}

setMethod("dim", signature("mixCopula"), function(x) dim(x@cops[[1]]))

## get parameters
setMethod("getParam", signature("mixCopula", "logicalOrMissing", "logicalOrMissing",
				"logicalOrMissing"),
	  function(copula, freeOnly = TRUE, attr = FALSE, named = attr) {
	      d <- dim(copula)
	      if(attr) { # more structured result
		  parC <- lapply(copula@cops, getParam,
				 freeOnly=freeOnly, attr=TRUE, named=named)
		  ## need attr(*, "param.(up|low)bnd") for for loglikCopula()
		  lowb <- lapply(parC, attr, "param.lowbnd")
		  uppb <- lapply(parC, attr, "param.upbnd")
		  structure(unlist(parC),
			    param.lowbnd = unlist(lowb),
			    param.upbnd  = unlist(uppb))
	      } else {
		  parC <- unlist(lapply(copula@cops, getParam,
					freeOnly=freeOnly, attr=FALSE, named=named))
		  c(parC, copula@w) ## <<- FIXME re-parametrize a la nor1mix:: (??)
	      }
	  })

setMethod("freeParam<-", signature("mixCopula", "numeric"),
	  function(copula, value) {
	      cops <- copula@cops
	      m <- length(cops)
	      nj <- vapply(cops, nFreeParam, 1)
	      ## FIXME re-parametrize a la nor1mix::nM2par / .par2nM
	      ## ----- i.e. value would only contain  qlogis(w[-1])  !!
	      nw <- length(iF.w <- isFree(w <- copula@w))
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

setMethod("pCopula", signature("numeric", "mixCopula"), pMixCopula)
setMethod("pCopula", signature("matrix",  "mixCopula"), pMixCopula)

##' The  c() function :
dMixCopula <- function(u, copula, log = FALSE, ...) {
    as.vector(
	vapply(copula@cops, dCopula, FUN.VALUE=numeric(nrow(u)), u=u, log=log, ...)
	%*%
	copula@w)
}

setMethod("dCopula", signature("numeric", "mixCopula"), dMixCopula)
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
