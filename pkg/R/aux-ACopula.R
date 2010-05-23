#### Functions and Methods for "ACopula" objects
#### class definition in ./AllClass.R

### ==== Clayton ===============================================================

##' Note: this is mainly to show that this function can be very well
##' approximated much more simply by just using m <- round(V0).
##' @param V0 numeric vector >= 0
##' @return optimal constant m for the fast rejection algorithm
##' @author Martin Maechler (based on Marius Hofert's code)
m.opt.retst <- function(V0){
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

##' Sample a random variate St ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)
##' *V_0)^{1/alpha}, V_0*I_{alpha = 1}, I_{alpha != 1}; 1) with
##' Laplace-Stieltjes transform exp(-V_0((1+t)^alpha-1)), see Nolan's book for
##' the parametrization, via an m-fold sum of random variates from
##' \tilde{S}(alpha, 1, (cos(alpha*pi/2)*V_0/m)^{1/alpha}, (V_0/m)
##' *I_{alpha = 1}, I_{alpha != 1}; 1) with Laplace-Stieltjes transform
##' exp(-(V_0/m)*((1+t)^alpha-1)).
##' @param m number of summands, any positive integer
##' @param V0 random variate
##' @param alpha parameter in (0,1]
##' @return St
##' @author Marius Hofert, Martin Maechler
retstablerej <- function(m,V0,alpha){
    gamm. <- (cos(alpha*pi/2)*V0/m)^(1/alpha)
    sum(unlist(lapply(integer(m),
                      function(.) {
                          ## apply standard rejection for sampling
                          ## \tilde{S}(alpha, 1, (cos(alpha*pi/2)
                          ##	*V_0/m)^{1/alpha}, (V_0/m)*I_{alpha = 1},
                          ## h*I_{alpha != 1}; 1) with Laplace-Stieltjes
                          ## transform exp(-(V_0/m)*((h+t)^alpha-h^alpha))
                          repeat {
                              V__ <- rstable1(1, alpha, beta=1, gamma = gamm.)
                              if(runif(1) <= exp(-V__))
                                  return(V__)
                          }})
               ## on acceptance, St_k ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)
               ## *V_0/m)^{1/alpha}, (V_0/m)*I_{alpha = 1}, h*I_{alpha != 1};
               ## 1) with Laplace-Stieltjes transform
               ## exp(-(V_0/m)*((h+t)^alpha-h^alpha))
               ))
}

##' Sample a vector of random variates St ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)
##' *V_0)^{1/alpha}, V_0*I_{alpha = 1}, I_{alpha != 1}; 1) with
##' Laplace-Stieltjes transform exp(-V_0((1+t)^alpha-1)), see Nolan's book for
##' the parametrization. This procedure calls retstablerej.
##' @param alpha parameter in (0,1]
##' @param V0 vector of random variates
##' @param h non-negative real number
##' @return vector of variates St
##' @author Marius Hofert, Martin Maechler
retstableR <- function(alpha,V0, h = 1){
    n <- length(V0)
    stopifnot(n >= 1, is.numeric(alpha), length(alpha) == 1,
              0 <= alpha, alpha <= 1) ## <- alpha > 1 ==> cos(pi/2 *alpha) < 0
    ## case alpha==1
    if(alpha == 1) { # alpha == 1 => St corresponds to a point mass at V0 with
        return(V0) # Laplace-Stieltjes transform exp(-V0*t)
    }
    ## else alpha != 1 : call fast rejection algorithm with optimal m
    m <- m.opt.retst(V0)
    mapply(retstablerej, m=m, V0=V0, alpha=alpha)
}

##' Sample a vector of random variates St ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)
##' *V_0)^{1/alpha}, V_0*I_{alpha = 1}, I_{alpha != 1}; 1) with
##' Laplace-Stieltjes transform exp(-V_0((1+t)^alpha-1)), see Nolan's book for
##' the parametrization. This procedure is more efficient than retstableR since
##' it calls the C function restable_c.
##' @param alpha parameter in (0,1]
##' @param V0 vector of random variates
##' @param h non-negative real number
##' @param method which method to call ("Marius Hofert", "Luc Devroye")
##' @return vector of variates St
##' @author Martin Maechler
## FIXME: find sophisticated default-method
retstableC <- function(alpha, V0, h = 1, method = c("MH","LD")){
    n <- length(V0)
    stopifnot(n >= 1, is.numeric(alpha), length(alpha) == 1,
              0 < alpha, alpha <= 1,
              is.numeric(h), length(h) == 1, h > 0)
    if(alpha == 1) { # alpha == 1 => St corresponds to a point mass at V0 with
	V0 # Laplace-Stieltjes transform exp(-V0*t)
    }
    else {
	method <- match.arg(method)
	.Call(retstable_c, V0, h = h, alpha, method)
    }
}

## switch to make fast C code the default
if(FALSE)
    retstable <- retstableR
retstable <- retstableC

### ==== Frank =================================================================

##' Random number generator for a Log(p) distribution, R version
##' @param n number of random variates to be generated
##' @param p parameter in (0,1)
##' @return vector of random variates from Log(p)
##' @author Marius Hofert, Martin Maechler
rlogR <- function(n,p) {
    stopifnot((n <- as.integer(n)) >= 0, 0 < p, p < 1)
    vec <- numeric(n)
    if(n >= 1) {
        u <- runif(n)
        l1 <- u > p
        vec[l1] <- 1
        i2 <- which( !l1 ) # of shorter length, say n2
        q2 <- 1-(1-p)^runif(length(i2)) # length n2
        l3 <- u[i2] < q2*q2
        i3 <- i2[l3]
        vec[i3] <- floor(1+abs(log(u[i3])/log(q2[l3])))
        l4 <- u[i2] > q2
        vec[i2[l4]] <- 1
        l5 <- ! (l3 | l4) # q2^2 <= u[i2] <= q2
        vec[i2[l5]] <- 2
    }
    vec
}

##' Random number generator for a Log(p) distribution, C version
##' @param n number of random variates to be generated
##' @param p parameter in (0,1)
##' @return vector of random variates from Log(p)
##' @author Martin Maechler
rlog <- function(n,p) {
    stopifnot(n >= 0,  0 < p, p < 1)
    .Call(rLog_c, n, p)
}

##' Sample V ~ F with Laplace-Stieltjes transform
##' (1-(1-exp(-t)*(1-e^(-theta1)))^alpha)/(1-e^(-theta0))
##' via the algorithm of Hofert (2010). R version.
##' @param p parameter 1-e^(-theta1)
##' @param alpha parameter theta0/theta1 in (0,1]
##' @param theta0_le_1 in {0,1} with 1 if and only if theta0 <= 1
##' @return V
##' @author Marius Hofert, Martin Maechler
rejFFrankR <- function(p,alpha,theta0_le_1) {
    if(theta0_le_1) {
	repeat{
	    U <- runif(1)
	    X <- rlog(1,p)
	    if(U*(X-alpha) <= 1/beta(X,1-alpha)) break
	}
    } else {
	repeat{
	    U <- runif(1)
	    X <- rFJoe(1,alpha)
	    if(U <= p^(X-1)) break
	}
    }
    X
}

##' Vectorize rejFFrankR. Generate a vector of variates V ~ F with
##' Laplace-Stieltjes transform (1-(1-exp(-t)*(1-e^(-theta1)))^alpha)
##' /(1-e^(-theta0)). R version.
##' @param n length of the vector of random variates
##' @param theta0 parameter theta0 in (0,infinity)
##' @param theta1 parameter theta1 in [theta0, infinity)
##' @return vector of random variates from F
##' @author Martin Maechler
rFFrankR <- function(n,theta0,theta1) {
    sapply(rep.int(-expm1(-theta1), n), rejFFrankR,
           alpha = theta0/theta1,
           theta0_le_1 = (theta0 <= 1))
}

##' Generate a vector of variates V ~ F with Laplace-Stieltjes transform
##' (1-(1-exp(-t)*(1-e^(-theta1)))^alpha)/(1-e^(-theta0)). C version.
##' @param n length of the vector of random variates
##' @param theta0 parameter theta0 in (0,infinity)
##' @param theta1 parameter theta1 in [theta0, infinity)
##' @return vector of random variates from F
##' @author Martin Maechler
rFFrank <- function(n,theta0,theta1) {
    .Call(rFFrank_c, n, theta0, theta1);
}

### ==== Joe ===================================================================

##' Sample V ~ F with F(n) = 1-1/(n*B(n,1-alpha)), n in IN, with
##' Laplace-Stieltjes transform 1-(1-exp(-t))^alpha via the algorithm of
##' Hofert (2010). R version.
##' @param n  sample size
##' @param alpha parameter
##' @return vector of random variates V
##' @author Marius Hofert, Martin Maechler
rFJoeR <- function(n,alpha) {
    stopifnot((n <- as.integer(n)) >= 0)
    V <- numeric(n)
    if(n >= 1) {
        if(alpha == 1) {
            V[] <- 1
        } else {
            u <- runif(n)
            ## FIXME(MM): (for alpha not too close to 1): re-express using 1-u
            l1 <- u <= alpha
            V[l1] <- 1
            i2 <- which(!l1)
            Ginv <- ((1-u[i2])*gamma(1-alpha))^(-1/alpha)
            floorGinv <- floor(Ginv)
            l3 <- (1-1/(floorGinv*beta(floorGinv,1-alpha)) < u[i2])
            V[i2[l3]] <- ceiling(Ginv[l3])
            i4 <- which(!l3)
            V[i2[i4]] <- floorGinv[i4]
        }
    }
    V
}

##' Sample V ~ F with F(n) = 1-1/(n*B(n,1-alpha)), n in IN, with
##' Laplace-Stieltjes transform 1-(1-exp(-t))^alpha via the algorithm of
##' Hofert (2010). C version.
##' @param n  sample size
##' @param alpha parameter
##' @return vector of random variates V
##' @author Martin Maechler
rFJoe <- function(n,alpha) {
    stopifnot(is.numeric(n), n >= 0)
    .Call(rFJoe_c, n, alpha)
}

### ==== other stuff ===========================================================

##' Function for setting the parameter in an ACopula
##' @param x ACopula
##' @param value parameter value
##' @return ACopula with theta set to value
##' @author Martin Maechler
setTheta <- function(x, value){
    stopifnot(is(x, "ACopula"),
	      is.na(value) || is.numeric(value))
    if(x@paraConstr(value)) ## parameter constraints are fulfilled
	x@theta <- value
    else
	stop("theta (=", format(value), ")  does not fulfil paraConstr()")
    x
}

## such that a user can have, e.g., myCop <- copAMH; theta(myCop) <- 3
if(FALSE) { ## this is just "didactical", not really useful
    setGeneric("theta<-", function(x, value) standardGeneric("theta<-"))
    setMethod("theta<-", "ACopula", setTheta)
}

##' Construct "paraConstr" function from an "interval"
##' @param int interval
##' @return parameter constraint function
##' @author Martin Maechler
mkParaConstr <- function(int){
    stopifnot(is(int, "interval")) # for now
    is.o <- int@open
    eL <- substitute(LL <= theta, list(LL = int[1])); if(is.o[1]) eL[[1]] <-
	as.symbol("<")
    eR <- substitute(theta <= RR, list(RR = int[2])); if(is.o[2]) eR[[1]] <-
	as.symbol("<")
    bod <- substitute(length(theta) == 1 && LEFT && RIGHT,
                      list(LEFT = eL, RIGHT= eR))
    as.function(c(alist(theta=), bod), parent.env(environment()))
    ## which is a fast version of
    ## r <- function(theta) {}
    ## environment(r) <- parent.env(environment())
    ## body(r) <- bod
    ## r
}
