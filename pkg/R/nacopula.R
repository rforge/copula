#### Implementation of function evaluations and random number generation
#### for nested Archimedean copulas

### returns the copula value at a certain vector u 
#FIXME: maybe make this applicable to a matrix of u's?
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

### returns U : matrix(*,n,d)
setGeneric("rn", function(x,n) standardGeneric("rn"))

setMethod("rn", signature(x ="final_nACopula"),
    function(x,n) {
	      Cout <- x@copula#outer copula
	      theta <- Cout@theta#theta for outer copula
	      V0 <- Cout@V0(n,theta)#generate V0's
	      childdat <- lapply(x@childCops,rnchild)#start recursion
	      dns <- length(x@comp)#dimension of the non-sectorial part
        r <- matrix(runif(n*dns),nrow=n,ncol=dns)#generate the non-sectorial part
        mat <- cbind(r, do.call(cbind, childdat$U))#put pieces together
        dat <- Cout@psi(-log(mat)/V0,theta=theta)#transform
        j <- c(x@comp,childdat$indCol)#get correct sorting order
        Cout@psi(dat[,j],theta)#permute data and return
    })

### returns list(U = matrix(*,n,d), indCol = vector of length d)
setGeneric("rnchild", function(x,n,psi0Inv,theta0,V0) standardGeneric("rnchild"))

### for all inner copulas
setMethod("rnchild", signature(x ="nACopula"),
    function(x,n,psi0Inv,theta0,V0) {
        Cin <- x@copula#inner copula
	      theta1 <- Cin@theta#theta_1 for inner copula
	      V01 <- Cin@V01(n,theta0,theta1)#generate V01's
	      childdat <- lapply(x@childCops,rnchild,n=n,psi0Inv=Cin@psiInv,theta0=theta1,V0=V01)#recursion
	      dns <- length(x@comp)#dimension of the non-sectorial part
        r <- matrix(runif(n*dns),nrow=n,ncol=dns)#generate the non-sectorial part
        mat <- cbind(r, do.call(cbind, childdat$U))#put pieces together
        dat <- exp(-V*psi0Inv(C@psi(-log(mat),theta1),theta0))#transform
        j <- c(x@comp,childdat$indCol)#get correct sorting order
        list(U = dat, indCol = j)#get list and return
    })