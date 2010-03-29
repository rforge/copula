#### Implementation of function evaluations and random number generation
#### for nested Archimedean copulas

### returns the copula value at a certain vector u
#FIXME: maybe make this applicable to a matrix of u's?
setGeneric("value", function(x, u) standardGeneric("value"))

setMethod("value", signature(x ="nACopula"),
    function(x,u) {
	stopifnot(is.numeric(u),
		  length(u) >= dim(x))	# can be larger
	C <- x@copula
	th <- C@theta
	## Now use u[j] for the direct components 'comp'
	C@psi(sum(unlist(lapply(u[x@comp], C@psiInv, theta=th)),
		  C@psiInv(unlist(lapply(x@childCops, value, u = u)),
			   theta=th)),
	      theta=th)
    })

### returns U : matrix(*,n,d)
setGeneric("rn", function(x,n,...) standardGeneric("rn"))

setMethod("rn", signature(x = "outer_nACopula"),
	  function(x, n, ...)
      {
	  Cp <- x@copula		# outer copula
	  theta <- Cp@theta		# theta for outer copula
	  V0 <- Cp@V0(n,theta)		# generate V0's
	  childL <- lapply(x@childCops, rnchild, # <-- start recursion
			     n=n,psi0Inv=Cp@psiInv,theta0=theta,V0=V0,...)
	  dns <- length(x@comp)	 # dimension of the non-sectorial part
	  r <- matrix(runif(n*dns), n, dns) # generate the non-sectorial part
	  mat <- cbind(r, do.call(cbind,lapply(childL, `[[`, "U"))) # put pieces together
	  mat <- Cp@psi(-log(mat)/V0, theta=theta) # transform
          ## get correct sorting order:
	  j <- c(x@comp, unlist(lapply(childL, `[[`, "indCol")))
	  ## extra checks -- comment later:
	  stopifnot(length(j) == ncol(mat))
	  m <- Cp@psi(mat[,j], theta)	  # permute data and return
	  ## extra checks:
	  stopifnot(length(dm <- dim(m)) == 2, dm == dim(mat))
	  m
      })

### returns list(U = matrix(*,n,d), indCol = vector of length d)
setGeneric("rnchild", function(x, n, psi0Inv, theta0, V0, ...) standardGeneric("rnchild"))

### for all inner copulas
setMethod("rnchild", signature(x ="nACopula"),
	  function(x, n, psi0Inv, theta0, V0,...)
      {
	  Cp <- x@copula		# inner copula
	  ## Consistency checks -- for now {comment later} :
	  stopifnot(is(Cp, "ACopula"), is.numeric(n), n == as.integer(n),
		    is.function(psi0Inv), is.numeric(V0), length(V0) == n,
		    is.numeric(theta0))
	  theta1 <- Cp@theta		    # theta_1 for inner copula
	  V01 <- Cp@V01(V0, theta0,theta1,...)	 # generate V01's
	  childL <- lapply(x@childCops, rnchild, # <-- recursion
			     n=n, psi0Inv = Cp@psiInv, theta0=theta1, V0=V01,...)
	  dns <- length(x@comp)	 # dimension of the non-sectorial part
	  r <- matrix(runif(n*dns), n, dns) # generate the non-sectorial part
	  ## put pieces together: first own comp.s, then the children's :
	  mat <- cbind(r, do.call(cbind, lapply(childL, `[[`, "U")))
	  mat <- exp(-V0* psi0Inv(Cp@psi(-log(mat), theta1),
				  theta0)) # transform
          ## get correct sorting order:
	  j <- c(x@comp, unlist(lapply(childL, `[[`, "indCol")))
	  list(U = mat, indCol = j)	   # get list and return
      })

if(FALSE) { ## evaluate the following into your R session if you need debugging:
trace(rn,      browser, exit=browser, signature=signature(x ="outer_nACopula"))

trace(rnchild, browser, exit=browser, signature=signature(x ="nACopula"))

}
