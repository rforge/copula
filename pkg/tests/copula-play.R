library(nacopula)
set.seed(1)

if(!dev.interactive())
    pdf("copula-play.pdf")

##====testing psi====

myCop <- setTheta(copAMH, value = 0.5) # is maybe more natural

setGeneric("psi", function(cop) standardGeneric("psi"))
setMethod(psi, "ACopula",
          function(cop) { function(t) cop@psi(t, theta = cop@theta) })
psi(myCop) ## is a *function*
psi(myCop)(0:4)
curve(psi(myCop)(x), 0, 4)
## but this can also be done directly [ => same curve "on top" :]
curve(myCop@psi(x, theta = myCop@theta),  0, 4, col = 2, add = TRUE)

##====testing Kendall's tau====

p.Tau <- function(cop, n = 201, xlim = pmin(paraI, 50), ...) {
    stopifnot(is(cop, "ACopula"))
    paraI <- cop@paraInterval
    theta <- seq(xlim[1], xlim[2], length.out = n)
    tit <- substitute(tau[NAME](theta), list(NAME = cop@name))
    plot(theta, cop@tau(theta), type = "l", main = tit, ...)
}

p.Tau(copAMH)
p.Tau(copClayton)
p.Tau(copFrank)
p.Tau(copGumbel)
p.Tau(copJoe)

##====test function====

tstCop <- function(cop, theta1 = cop@theta,
                   thetavec = cop@theta,
                   i10 = 1:10, nRnd = 50,
                   t01 = (1:63)/64, ## exact binary fractions
                   lTDCvec = NA_real_, uTDCvec = NA_real_)
{
    stopifnot(is(cop, "ACopula"))
    cat0 <- function(...) cat(..., "\n", sep = "")
    theta0 <- cop@theta
    CT <- list()
    cat0(sprintf("(1) copula family: %10s, theta0 = %g",
                 cop@name, theta0))
    cat("\n(2) values of psi at i10:\n")
    CT <- c(CT, list(psi = system.time(p.i <- cop@psi(i10,theta = theta0))))
    print(p.i)
    cat("check if psi(Inf)=0: ")
    stopifnot(cop@psi(Inf, theta = theta0)==0)
    cat0("TRUE")
    cat("check if psiInv(numeric(0)) is numeric(0): ")
    n0 <- numeric(0)
    stopifnot(identical(n0, cop@psiInv(n0, theta = theta0)))
    cat0("TRUE")
    cat("check if psiInv(0)=Inf: ")
    stopifnot(cop@psiInv(0, theta = theta0)==Inf)
    cat0("TRUE")
    cat0("values of psiInv at t01:\n")
    CT <- c(CT, list(psiI = system.time(pi.t <- cop@psiInv(t01,theta = theta0))))
    print(pi.t)
    CT[["psiI"]] <- CT[["psiI"]] + system.time(pi.pi <- cop@psiInv(p.i, theta = theta0))
    CT[["psi" ]] <- CT[["psi" ]] + system.time(p.pit <- cop@psi(pi.t,  theta = theta0))
    cat0("check if psiInv(psi(i10))==i10: ", all.equal(pi.pi, i10))
    cat0("check if psi(psiInv(t01))==t01: ", all.equal(p.pit, t01))
    cat("\n(3) parameter interval:\n")
    print(cop@paraInterval)
    cat0("theta1=",theta1)
    cat0("nesting condition for theta0 and theta1 fulfilled: ",
         cop@nestConstr(theta0,theta1))
    CT <- c(CT, list(V0 = system.time(V0 <- cop@V0(nRnd,theta0))))
    cat0("\n(4) ",nRnd," generated V0's:")
    print(summary(V0))
    CT <- c(CT, list(V01 = system.time(V01 <- cop@V01(V0,theta0,theta1))))
    cat0(nRnd," generated V01's:")
    print(summary(V01))
    nt <- length(thetavec)
    cat("\n(5) tau at thetavec:\n")
    CT <- c(CT, list(tau = system.time(ta <- cop@tau(thetavec))))
    print(ta)
    CT <- c(CT, list(tauI = system.time(ta.I <- cop@tauInv(ta))))
    cat0("check if tauInv(tau(thetavec))==thetavec: ", all.equal(ta.I, thetavec))
    lTDCvec <- rep(as.double(lTDCvec), length.out= nt)
    uTDCvec <- rep(as.double(uTDCvec), length.out= nt)
    cat("\n(6) lTDC at thetavec:\n")
    CT <- c(CT, list(lTDC = system.time(lT <- cop@lTDC(thetavec))))
    CT <- c(CT, list(lT.I = system.time(lT.I <- cop@lTDCInv(lT))))
    print(lT)
    cat0("check if lTDCInv(lTDC(thetavec))==lTDCvec: ", all.equal(lT.I, lTDCvec))
    cat("\n(7) uTDC at thetavec:\n")
    CT <- c(CT, list(uTDC = system.time(uT <- cop@uTDC(thetavec))))
    CT <- c(CT, list(uT.I = system.time(uT.I <- cop@uTDCInv(uT))))
    print(uT)
    cat0("check if uTDCInv(uTDC(thetavec))==uTDCvec: ", all.equal(uT.I, uTDCvec))
    ## return the list of CPU-time measurements
    ## as a classed object with a convenient print method:
    class(CT) <- "proc_time_list"
    CT
}

## print() method for the tstCop() results
print.proc_time_list <- function (x, ...) {
    stopifnot(is.list(x), !is.null(nx <- names(x)))
    for(nm in nx)
	if(!all(x[[nm]] == 0)) {
	    cat(nm,":\n"); print(x[[nm]], ...)
	}
    invisible(x)
}

##====copAMH====

## define object
myAMH <- setTheta(copAMH, 0.7135001)
thetavec <- c(0.1,0.3,0.5,0.7,0.9)

## call test function
tstCop(myAMH, 0.9429679, thetavec = thetavec)

##====copClayton====

## define object
myAMH <- setTheta(copAMH, 0.7135001)

myClayton <- setTheta(copClayton, 0.5)
thetavec <- c(0.5,1,2,5,10)

## call test function
tstCop(myClayton, 2, thetavec, lTDC = thetavec, uTDC = NA)

##====copFrank===

## define object
myFrank <- setTheta(copFrank, 1.860884)
thetavec <- c(0.5,1,2,5,10)

## call test function
tstCop(myFrank, 5.736283, thetavec)

##====copGumbel===

## define object
myGumbel <- setTheta(copGumbel, 1.25)
thetavec <- c(1,2,4,6,10)

## call test function
tstCop(myGumbel,2, thetavec, lTDC = NA, uTDC = thetavec)

##====copJoe===

## define object
myJoe <- setTheta(copJoe, 1.25)
thetavec <- c(1.1,2,4,6,10)

## call test function
tstCop(myJoe, 2, thetavec, lTDC = NA, uTDC = thetavec)
