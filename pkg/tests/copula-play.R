library(nacopula)

if(!dev.interactive())
    pdf("copula-play.pdf")

myCop <- copAMH
theta(myCop) <- 0.5

cop2 <- setTheta(copAMH, value = 0.5) # is maybe more natural
stopifnot( identical(cop2, myCop) )

setGeneric("psi", function(cop) standardGeneric("psi"))
setMethod(psi, "ACopula",
          function(cop) { function(t) cop@psi(t, theta = cop@theta) })
psi(myCop) ## is a *function*
psi(myCop)(0:4)
curve(psi(myCop)(x), 0, 4)
## but this can also be done directly [ => same curve "on top" :]
curve(myCop@psi(x, theta = myCop@theta),  0, 4, col = 2, add = TRUE)

## Kendall's tau --
(ii <- copClayton@paraInterval) # upper bound Inf
th <- seq(ii[1], 50, length = 201)
plot(th, copClayton@tau(th), type = "l",
     main = expression(tau["Clayton"](theta)))

p.Tau <- function(cop, n = 201, xlim = pmin(paraI, 50), ...) {
    stopifnot(is(cop, "ACopula"))
    paraI <- cop@paraInterval
    theta <- seq(xlim[1], xlim[2], length.out = n)
    tit <- substitute(tau[NAME](theta), list(NAME = cop@name))
    plot(theta, cop@tau(theta), type = "l", main = tit, ...)
}

p.Tau(copAMH)
p.Tau(copClayton)
if(FALSE) ## FIXME
p.Tau(copFrank)# -> error in integrate !
p.Tau(copGumbel)
p.Tau(copJoe)# << now works (somewhat slowly)

##====test function====

tstCop <- function(cop, theta1, thetavec,
                   lTDCvec = NA_real_, uTDCvec = NA_real_)
{
    stopifnot(is(cop, "ACopula"))
    i10 <- 1:10
    t01 <- (1:15)/16 ## exact binary fractions
    n0 <- numeric(0)
    n <- 20

    cat0 <- function(...) cat(..., "\n", sep = "")

    theta0 <- cop@theta
    cat0(sprintf("(1) copula family: %10s, theta0 = %g",
                 cop@name, theta0))
    cat("\n(2) values of psi at i10:\n")
    print(cop@psi(i10,theta = theta0))
    cat("psiInv(numeric(0)) must be n..(0):\n")
    stopifnot(identical(n0, cop@psiInv(n0, theta = theta0)))
    cat("values of psiInv at t01:\n")
    print(cop@psiInv(t01,theta = theta0))
    cat0("check if psiInv(psi(i10))==i10: ",
         all.equal(cop@psiInv(cop@psi(i10,theta = theta0),theta = theta0),i10))
    cat0("check if psi(psiInv(t01))==t01: ",
         all.equal(cop@psi(cop@psiInv(t01,theta = theta0),theta = theta0),t01))
    cat("\n(3) parameter interval:\n")
    print(cop@paraInterval)
    cat0("nesting condition for theta0 and theta1 fulfilled: ",
         cop@nestConstr(theta0,theta1))
    V0 <- cop@V0(n,theta0)
    cat0("\n(4) ",n," generated V0's:")
    print(V0)
    ## V01=cop@V01(V0,theta0,theta1)#todo: too slow for frank and joe
    ## cat(n," generated V01's:\n",sep="")
    ## print(V01)
    nt <- length(thetavec)
    cat("\n(5) tau at thetavec:\n")
    print(cop@tau(thetavec))
    cat0("check if tauInv(tau(thetavec))==thetavec: ",
        all.equal(cop@tauInv(cop@tau(thetavec)),thetavec))
    lTDCvec <- rep(as.double(lTDCvec), length.out= nt)
    uTDCvec <- rep(as.double(uTDCvec), length.out= nt)
    cat("\n(6) lTDC at thetavec:\n")
    print(cop@lTDC(thetavec))
    cat0("check if lTDCInv(lTDC(thetavec))==lTDCvec: ",
         all.equal(cop@lTDCInv(cop@lTDC(thetavec)), lTDCvec))
    cat("\n(7) uTDC at thetavec:\n")
    print(cop@uTDC(thetavec))
    cat0("check if uTDCInv(uTDC(thetavec))==uTDCvec: ",
         all.equal(cop@uTDCInv(cop@uTDC(thetavec)), uTDCvec))
    invisible()
}

##====test setup====

set.seed(1)

##====copAMH====

##define object
myAMH <- copAMH
theta(myAMH) <- 0.7135001
thetavec <- c(0.1,0.3,0.5,0.7,0.9)

##call test function
tstCop(myAMH, 0.9429679, thetavec = thetavec)

##====copClayton====

##define object
myClayton <- copClayton
theta(myClayton) <- 0.5
thetavec <- c(0.5,1,2,5,10)

##call test function
tstCop(myClayton, 2, thetavec, lTDC = thetavec, uTDC = NA)

##====copFrank===

##define object
myFrank <- copFrank
theta(myFrank) <- 1.860884
thetavec <- c(0.5,1,2,5,10)

tau.th <- c(0.055417, 0.11002, 0.21389, 0.4567, 0.66578)
tau.F <- myFrank@tau(thetavec)
stopifnot(all.equal(tau.th, tau.F, tol = 0.0001),
	  all.equal(myFrank@tauInv(tau.F),
		    thetavec, tol = 5e-5))# tauInv() using uniroot()

tstCop(myFrank, 5.736283, thetavec)

##====copGumbel===

##define object
myGumbel <- copGumbel
theta(myGumbel) <- 1.25
thetavec <- c(1,2,4,6,10)

##call test function
tstCop(myGumbel,2, thetavec, lTDC = NA, uTDC = thetavec)

##====copJoe===

##define object
myJoe <- copJoe
theta(myJoe) <- 1.25
thetavec <- c(1.1,2,4,6,10)

##call test function
tstCop(myJoe, 2, thetavec, lTDC = NA, uTDC = thetavec)
