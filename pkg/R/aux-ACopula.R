#### Functions and Methods for	"ACopula" objects (class definition in ./AllClass.R )

### ====general functions====

## standard rejection algorithm for sampling
## S(alpha,1,(cos(alpha*pi/2)V0/m)^(1/alpha),0;1) with
## Laplace-Stieltjes transform exp(-(V0/m)*t^alpha) 
retstablerej1 <- function(m,V0,alpha,h){#gets one integer m and one real V0>0
  repeat{
      Vtilde <- rstable1(1, alpha, beta=1,
                         gamma = (cos(alpha*pi/2)*V0/m)^(1/alpha))
      u <- runif(1)
      if(u <= exp(-h*Vtilde)) {
          variate <- Vtilde
          break
      }
  }
  variate
}

## standard rejection algorithm for sampling the m-fold sum of variates
## from the distribution with Laplace-Stieltjes transform
## exp(-(V0/m)*((h+t)^alpha-h^alpha))
## calls retstablerej1 m times for given m and evaluates the sum of the results
retstablerej <- function(m,V0,alpha,h){#gets one integer m and one real V0>0
  sum(unlist(lapply(rep(m,m),retstablerej1,V0=V0,alpha=alpha,h=h)))
}

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
        V0 # sample from S(1,1,0,V0;1) with Laplace-Stieltjes transform exp(-V0*t)
    }
    variates <- numeric(n)
    ## case alpha!=1 and h==0
    if(h == 0) {
      rstable1(n, alpha, beta=1,
               gamma = (cos(alpha*pi/2)*V0)^(1/alpha))
               ## sample from S(alpha,1,(cos(alpha*pi/2)V0)^(1/alpha),0;1)
               ## with Laplace-Stieltjes transform exp(-V0*t^alpha)
    }
    ## case alpha!=1 and h!=0
    ## call fast rejection algorithm
    ## determine optimal constant m for the fast rejection algorithm
    logc <- V0*h^alpha
    floorlogc <- floor(logc)
    ceilinglogc <- ceiling(logc)
    c <- exp(logc)
    floorvalue <- floorlogc*c^(1/floorlogc)
    ceilingvalue <- ceilinglogc*c^(1/ceilinglogc)
    m=numeric(n)
    l1=(logc<=1)
    i1=(1:n)[l1]
    m[i1]=1
    l2=!l1
    i2=(1:n)[l2]
    l3=(logc[i2]>1 && floorvalue[i2]<=ceilingvalue[i2])
    i3=i2[l3]
    m[i3]=floorlogc[i3]
    l4=!l3
    i4=i2[l4]
    m[i4]=ceilinglogc[i4]
    ## call rejection with these m
    mapply(retstablerej,m=m,V0=V0,alpha=alpha,h=h)
}

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
