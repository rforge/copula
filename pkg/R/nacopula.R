#### Implementation of function evaluations and random number generation
#### for nested Archimedean copulas

### returns the copula value at a certain vector u
#FIXME: maybe make this applicable to a matrix of u's?
setGeneric("value", function(x, u) standardGeneric("value"))

setMethod("value", signature(x ="nACopula"),
    function(x,u) {
	stopifnot(is.numeric(u), all(u)>=0, all(u)<=1,
		  length(u) >= dim(x))	# can be larger
	C <- x@copula
	th <- C@theta
	## Now use u[j] for the direct components 'comp'
	C@psi(sum(unlist(lapply(u[x@comp], C@psiInv, theta=th)),
		  C@psiInv(unlist(lapply(x@childCops, value, u = u)),
			   theta=th)),
	      theta=th)
    })
    
### returns the probability that a random vector following the given copula
### falls in the hypercube with lower and upper corner "l" and "u", 
### respectively
setGeneric("prob", function(x, l, u) standardGeneric("prob"))

setMethod("prob", signature(x ="outer_nACopula"),
    function(x, l,u) {
        d <- dim(x)
        ## TODO: maybe allow  l & u to be  k x d matrices
        ##        --> return vector of probabilities of length k
	stopifnot(is.numeric(l), is.numeric(u),
		  length(u) == d, d == length(l),
                  0 <= l, l <= u, u <= 1)
        D <- 2^d
        m <- 0:(D - 1)
        ## digitsBase() from package 'sfsmisc' {slightly simplified} :
        ## Purpose: Use binary representation of 0:N
        ## Author: Martin Maechler, Date:  Wed Dec  4 14:10:27 1991
        II <- matrix(0, nrow = D, ncol = d)
        for (i in d:1L) {
            II[,i] <- m %% 2L + 1L
            if (i > 1) m <- m %/% 2L
        }
        ## Sign: the ("u","u",...,"u") case has +1; = c(2,2,...,2)
        Sign <- c(1,-1)[1L + (- rowSums(II)) %% 2]
        U <- array(cbind(l,u)[cbind(c(col(II)), c(II))], dim = dim(II))
        sum(Sign * apply(U, 1, value, x=x))
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
	  m <- mat[,order(j)]	  # permute data and return
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
	  theta1 <- Cp@theta            # theta_1 for inner copula
          ## generate V01's (only for one sector since the
          ## recursion in rn() takes care of all sectors):
	  V01 <- Cp@V01(V0, theta0,theta1,...)
	  childL <- lapply(x@childCops, rnchild, # <-- recursion
                           n=n, psi0Inv = Cp@psiInv, theta0=theta1, V0=V01,...)
	  dns <- length(x@comp)	 # dimension of the non-sectorial part
	  r <- matrix(runif(n*dns), n, dns) # generate the non-sectorial part
	  ## put pieces together: first own comp.s, then the children's :
	  mat <- cbind(r, do.call(cbind, lapply(childL, `[[`, "U")))
	  mat <- exp(-V0* psi0Inv(Cp@psi(-log(mat)/V01, theta1),
				  theta0)) # transform
          ## get correct sorting order:
	  j <- c(x@comp, unlist(lapply(childL, `[[`, "indCol")))
	  list(U = mat, indCol = j)     # get list and return
      })

if(FALSE) { ## evaluate the following into your R session if you need debugging:
trace(rn,      browser, exit=browser, signature=signature(x ="outer_nACopula"))

trace(rnchild, browser, exit=browser, signature=signature(x ="nACopula"))

}

###--- NACopula constructor

##' <description>
##'
##' <details>
##' @title THE outer_nACopula constructor function
##' @param family character string: short or longer form of copula family name
##' @param nACform a "formula" of the form C(th, comp, list(C(..), C(..)))
##' @return a valid outer_nACopula object
##' @author Martin Maechler
onACopula <- function(family, nACform) {
    nacl <- substitute(nACform)
    stopifnot(is.character(family), length(family) == 1,
              identical(nacl[[1]], as.symbol("C")))
    nacl[[1]] <- as.symbol("oC")
    if(nchar(family) <= 2)# it's a short name
        family <- c_longNames[family]
    stopifnot(is(COP <- get(c_objNames[family]),#, envir = "package:nacopula"
                 "ACopula"))
    mkC <- function(cClass, a,b,c) {
        if(missing(b) || length(b) == 0) b <- integer()
        if(missing(c) || length(c) == 0) c <- list()
        else if(length(c) == 1 && !is.list(c)) c <- list(c)
        else stopifnot(is.list(c))
        stopifnot(is.numeric(a), length(a) == 1, is.numeric(b))
        if(any(sapply(c, class) != "nACopula"))
            stop("third entry of 'nACform' must be NULL or) list of 'C(..)' terms")
        new(cClass, copula = setTheta(COP, a),
            comp = as.integer(b), childCops = c)
    }
    C <- function(a,b,c) mkC("nACopula", a,b,c)
    oC <- function(a,b,c) mkC("outer_nACopula", a,b,c)

    eval(nacl)
}
