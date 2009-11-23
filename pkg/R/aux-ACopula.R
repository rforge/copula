#### Functions and Methods for	"ACopula" objects (class definition in ./AllClass.R )

### ====general functions====

## compute random variates from an exponentially-tilted Stable distribution
## \tilde{S}(alpha,1,(cos(alpha*pi/2)V0)^(1/alpha),V0*Indicator(alpha==1),h*Indicator(alpha!=1);1)
## with corresponding Laplace-Stieltjes transform exp(-V0*((h+t)^alpha-h^alpha))
retstable <- function(alpha,V0,h) {
    n <- length(V0)
    stopifnot(n >= 1,
              length(alpha) == 1, alpha >= 0,
              length(h) == 1, h >= 0)
    ## case alpha==1
    if(alpha == 1) {
        return(V0) # sample from S(1,1,0,V0;1) with Laplace-Stieltjes transform exp(-V0*t)
    }
    variates <- numeric(n)
    ## case alpha!=1 and h==0
    if(h == 0) {
        for(i in 1:n) {
            ## sample from S(alpha,1,(cos(alpha*pi/2)V0)^(1/alpha),0;1)
            ## with Laplace-Stieltjes transform exp(-V0*t^alpha)
            variates[i] <- rstable1(1, alpha, beta=1,
                                    gamma = (cos(alpha*pi/2)*V0[i])^(1/alpha))
        }
        return(variates)
    }
    ## FIXME: move parts out of loop below:
    ## determine optimal constant m (*vector*) for the fast rejection algorithm

    ## case alpha!=1 and h!=0
    ## call fast rejection algorithm
    for(i in 1:n) {
        ## determine optimal constant m for the fast rejection algorithm
        logc <- V0[i]*h^alpha
        floorlogc <- floor(logc)
        ceilinglogc <- ceiling(logc)
        c <- exp(logc)
        floorvalue <- floorlogc*c^(1/floorlogc)
        ceilingvalue <- ceilinglogc*c^(1/ceilinglogc)
        if(logc <= 1) {
            m <- 1
        } else if(logc > 1 && floorvalue <= ceilingvalue) {
            m <- floorlogc
        } else {
            m <- ceilinglogc
        }
        ## apply standard rejection algorithm to sample the m-fold sum of variates
        ## from the distribution with Laplace-Stieltjes transform
        ## exp(-(V0/m)*((h+t)^alpha-h^alpha))
        variates[i] <- 0
        for(j in 1:m) {
            repeat{
                ## sample from S(alpha,1,(cos(alpha*pi/2)V0/m)^(1/alpha),0;1) with
                ## Laplace-Stieltjes transform exp(-(V0/m)*t^alpha)
                Vtilde <- rstable1(1, alpha, beta=1,
                                   gamma = (cos(alpha*pi/2)*V0[i]/m)^(1/alpha))
                u <- runif(1)
                if(u <= exp(-h*Vtilde)) {
                    variates[i] <- variates[i]+Vtilde
                    break
                }
            }
        }
    }
    return(variates)
}



## such that a useR can have	myCop <- copAMH; theta(myCop) <- 3  :

setGeneric("theta<-", function(x, value) standardGeneric("theta<-"))

setTheta <- function(x, value) {
    stopifnot(is(x, "ACopula"),
	      is.na(value) || is.numeric(value))
    if(x@paraConstr(value)) ## parameter constraints are fulfilled
	x@theta <- value
    else
	stop("theta (=", format(value), ")  does not fulfil paraConstr()")
    x
}


setMethod("theta<-", "ACopula", setTheta)


mkParaConstr <- function(int) {
    ## Purpose: Construct  'paraConstr' function from an "interval"
    ## --------------------------------------------------------------
    ## Author: Martin Maechler, Date: 16 Nov 2009, 12:40
    stopifnot(is(int, "interval"))	# for now
    is.o <- int@open
    eL <- substitute(LL <= theta, list(LL = int[1])); if(is.o[1]) eL[[1]] <- as.symbol("<")
    eR <- substitute(theta <= RR, list(RR = int[2])); if(is.o[2]) eR[[1]] <- as.symbol("<")
    bod <- substitute(length(theta) == 1 && LEFT && RIGHT,
                      list(LEFT = eL, RIGHT= eR))
    ##
    as.function(c(alist(theta=), bod), parent.env(environment()))
    ## which is a fast version of
    ##
    ## r <- function(theta) {}
    ## environment(r) <- parent.env(environment())
    ## body(r) <- bod
    ## r
}
