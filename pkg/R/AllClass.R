
###  ???? Should we try to be compatible to  'Copula' package ??

### Mother class of all (simple, *NON*-nested) Archimedean Copula Types
### for *any* dimension d
setClass("ACopula",
	 representation(name = "character",
                        psi = "function",    # of (t, theta) -- the generator
                        psiInv = "function", # of (p, theta) -- psi_inverse: \psi^{-1}(p) = t
                        paraConstr = "function", # of (theta) ; constr(theta) |--> TRUE: "fulfilled"
                        V0 = "function",	# of (n,theta) -- RNGenerator
                        tau= "function",	# of (theta)
                        tauInv = "function",    # of (tau)
                        lTDC = "function",    # of (theta) lower bound  \lambda_l
                        lTDCInv = "function", # of (lambda) - Inverse of \lambda_l
                        uTDC = "function",    # of (theta)  - upper bound  \lambda_u
                        uTDCInv = "function", # of (lambda) - Inverse of \lambda_u

                        ## Nesting properties if the child copulas are of the same family :
                        nestConstr = "function", # of (th0, th1) ; TRUE <==> "fulfilled"
                        V01= "function"	# of (V0,theta0,theta1)
                        ),
         ## prototype = prototype(Dim = integer(2), Dimnames = list(NULL,NULL)),
	 validity = function(object) {
	     if (length(nm <- object@name) != 1 || nchar(nm) < 1)
		 return("'name' must be a string (of at least 1 character)")
             checkFun <- function(sName, nArgs) {
                 f <- slot(object, sName)
                 if (length(formals(f)) == nArgs)
                     TRUE
                 else
                     paste("slot '",sName,"' must be a function of ",nArgs,
                           " arguments", sep="")
             }
             if(!isTRUE(tt <- checkFun("psi", 2)))	return(tt)
             if(!isTRUE(tt <- checkFun("psiInv", 2)))	return(tt)
             if(!isTRUE(tt <- checkFun("paraConstr", 1))) return(tt)
             if(!isTRUE(tt <- checkFun("nestConstr", 2))) return(tt)
             ## ...
             ## ... (TODO)

             ## Check more :
	     if (object@psi(0, theta= 1/2) != 1)
		 return("psi(0, th=1/2) != 1 -- seems an invalid generator")

             ## ....
             ## ....
             ## ....
             ## ....

	     ## 'else'	ok :
	     TRUE
	 })


### Nested Archimedean Copulas with *specified* dimension(s)
setClass("nACopula",
	 representation(copula = "ACopula",
                        comp = "integer", # from 1:d -- of length in [0,d]
                        childCops = "list" #of nACopulas, possibly empty
                        ## TODO? nesting properties (V01,..) for specific mother-child relations
                        ),
         validity = function(object) {
             if(length(d <- dim(object)) != 1 || !is.numeric(d) || d <= 0)
                 return("invalid dim(.)")
             ic <- object@comp
             if((lc <- length(ic)) > d)
                 return("more components than the dimension")
             if(any(ic < 1) || any(ic > d))
                 return(paste("components (indices) must be in 1:d, d=",d))
             if(!all(ok <- "nACopula" == sapply(object@childCops, class)))
                 return("All 'childCops' elements must be 'nACopula' objects")
             iChilds <- unlist(lapply(object@childCops, comp))
             allC <- sort(c(ic, iChilds))
             if(length(allC) != d)
                 return("must have d coordinates (from 'comps' and child copulas)")
             if(!all(allC == 1:d))
                 return(paste("The implicit coordinates are not identical to 1:d; instead\n  ",
                              paste(allC, collapse=", ")))
             ##
             TRUE
         })


## The dim() method is nicely defined  *recursive*ly :
setMethod("dim", signature(x = "nACopula"),
	  function(x) length(x@comp) + sum(sapply(x@childCops, dim)))

