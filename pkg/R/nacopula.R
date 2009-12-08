#### Implementation of function evaluations and random number generation
#### for nested Archimedean copulas

###
setGeneric("value", function(x, u) standardGeneric("value"))

setMethod("value", signature(x ="nACopula"),
    function(x,u) {
	stopifnot(is.numeric(u),
		  length(u) >= dim(x)) # can be larger
	C <- x@copula
	th <- C@theta
        ## Now use u[j] for the direct components 'comp'
	C@psi(sum(unlist(lapply(u[x@comp], C@psiInv, theta=th)),
		  C@psiInv(unlist(lapply(x@childCops, value, u = u)),
                           theta=th)),
	      theta=th)
    })
