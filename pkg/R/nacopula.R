#### Implementation of function evaluations and random number generation of supported nested Archimedean copulas

### 
setGeneric("value", function(x, u) standardGeneric("value"))

setMethod("value", signature(x ="nACopula"),
  function(x,u)	{
      x@psi(sum(unlist(lapply(x@comp,x@psiInv,theta=x@theta)),
                x@psiInv(unlist(lapply(x@childCops,value)),theta=x@theta)),
            theta=x@theta)
}	)
	

