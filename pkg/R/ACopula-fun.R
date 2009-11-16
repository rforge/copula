#### Functions and Methods for	"ACopula" objects (class definition in ./AllClass.R )

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
