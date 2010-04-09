#### Functions and Methods for	"ACopula" objects (class definition in ./AllClass.R )

### ====general functions====


## standard rejection algorithm for sampling the m-fold sum of variates
## from the distribution with Laplace-Stieltjes transform
## exp(-(V0/m)*((1+t)^alpha-1))
## calls sampler  m times for given m and evaluates the sum of the results
retstablerej <- function(m,V0,alpha){ #gets one integer m and one real V0>0
  gamm. <- (cos(alpha*pi/2)*V0/m)^(1/alpha)
  sum(unlist(lapply(integer(m),
                    function(.) {
                        repeat {
                            V__ <- rstable1(1, alpha, beta=1, gamma = gamm.)
                            if(runif(1) <= exp(-V__))
                                return(V__)
                        }})
             ))
}

##' <description>
##' --- Note: this is mainly to show that this function can be very
##'     ====  well approximated much more simply :
##'           namely, just using   m <- round(V0)   !!
##' <details>
##' @title determine optimal constant m for the fast rejection algorithm
##' @param V0 numeric vector >= 0
##' @return integer vector, typically >= V0
##' @author Martin Maechler (based on Marius' code)
m.opt.retst <- function(V0)
{
    n <- length(V0)
    fV <- floor(V0)
    cV <- ceiling(V0)
    v1 <- fV*exp(V0/fV)
    v2 <- cV*exp(V0/cV)

    m <- integer(n)
    l1 <- (V0 <= 1)
    m[which(l1)] <- 1L

    i2 <- which(!l1) ## those with V0 > 1
    l3 <- (v1[i2] <= v2[i2])
    i3 <- i2[l3]
    m[i3] <- fV[i3]
    i4 <- i2[!l3]
    m[i4] <- cV[i4]
    m
}

## compute random variates from an exponentially-tilted Stable distribution
## \tilde{S}(alpha,1,(cos(alpha*pi/2)V0)^(1/alpha),V0*Indicator(alpha==1),Indicator(alpha!=1);1)
## with corresponding Laplace-Stieltjes transform exp(-V0*((1+t)^alpha-1))
retstableR <- function(alpha,V0) {
    n <- length(V0)
    stopifnot(n >= 1, is.numeric(alpha), length(alpha) == 1,
              0 <= alpha, alpha <= 1) ## <- alpha > 1 ==> cos(pi/2 *alpha) < 0
    ## case alpha==1
    if(alpha == 1) {
       return(V0) # sample from S(1,1,0,V0;1) with Laplace-Stieltjes transform exp(-V0*t)
    }
    ## else alpha != 1 : call fast rejection algorithm

    m <- m.opt.retst(V0)
    ## call rejection with these m
    mapply(retstablerej, m=m, V0=V0, alpha=alpha)
}

retstableC <- function(alpha, V0) {
    n <- length(V0)
    stopifnot(n >= 1, is.numeric(alpha), length(alpha) == 1,
              0 <= alpha, alpha <= 1) ## <- alpha > 1 ==> cos(pi/2 *alpha) < 0
    if(alpha == 1)
	V0 # sample from S(1,1,0,V0;1) with Laplace-Stieltjes transform exp(-V0*t)
    else
	.Call(retstable_c, V0, alpha)
    ## REAL PROBLEM: This gives *different* result than the pure R version !
}

## For now --- FIXME
retstable <- retstableR


setTheta <- function(x, value) {
    stopifnot(is(x, "ACopula"),
	      is.na(value) || is.numeric(value))
    if(x@paraConstr(value)) ## parameter constraints are fulfilled
	x@theta <- value
    else
	stop("theta (=", format(value), ")  does not fulfil paraConstr()")
    x
}

## such that a useR can have	myCop <- copAMH; theta(myCop) <- 3  :
if(FALSE) {## this is just "didactical", not really useful:
setGeneric("theta<-", function(x, value) standardGeneric("theta<-"))
setMethod("theta<-", "ACopula", setTheta)
}

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
