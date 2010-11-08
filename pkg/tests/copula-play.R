## Copyright (C) 2010 Marius Hofert and Martin Maechler
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

require(nacopula)
set.seed(1)

if(!dev.interactive())
    pdf("copula-play.pdf")

##' ==== testing psi ====

myCop <- setTheta(copAMH, value = 0.5) # is maybe more natural

setGeneric("psi", function(cop) standardGeneric("psi"))
setMethod(psi, "acopula",
          function(cop) { function(t) cop@psi(t, theta = cop@theta) })
psi(myCop) # is a function
psi(myCop)(0:4)
curve(psi(myCop)(x), 0, 4)
##' but this can also be done directly [ => same curve "on top" :]
curve(myCop@psi(x, theta = myCop@theta),  0, 4, col = 2, add = TRUE)

##' ==== testing Kendall's tau ====

p.Tau <- function(cop, n = 201, xlim = pmin(paraI, 50), ...) {
    stopifnot(is(cop, "acopula"))
    paraI <- cop@paraInterval
    theta <- seq(xlim[1], xlim[2], length.out = n)
    tit <- substitute(tau[NAME](theta), list(NAME = cop@name))
    plot(theta, cop@tau(theta), type = "l", main = tit, ...)
    abline(h = c(0,1), lty = 3, col = "gray20")
}

p.Tau(copAMH)
p.Tau(copClayton)
p.Tau(copFrank, xlim = c(0, 80), ylim= 0:1) # fast via debye_1()
p.Tau(copGumbel)
p.Tau(copJoe, ylim = 0:1, yaxs="i")

### ==== test function ==== ------------------------------------------------

##' @title stopifnot() plus output
##' @param expr
##' @param prefix
##' @param true
##' @return
##' @author Martin Maechler
checkifnot <- function(expr, prefix = "check if", true = "[Ok]")
{
    c0  <- function(...) cat(..., sep = "")
    ## match.call(): not "calling" expr too early:
    c0(prefix, deparse(match.call()[[2]])[1],": ")
    stopifnot(expr)
    c0(true,"\n")
}


##' @title Perform a set of checks on a copula object (with theta set)
##' @param cop copula
##' @param theta1 parameter theta1
##' @param thetavec vector of parameters
##' @param i10 values where psi is evaluated
##' @param nRnd number of generated V0's and V01's
##' @param t01 values where psiinv is evaluated
##' @param lambdaLvec vector of lower tail-dependence coefficients
##' @param lambdaUvec vector of upper tail-dependence coefficients
##' @return list of measurements
##' @author Marius Hofert, Martin Maechler
tstCop <- function(cop, theta1 = cop@theta, thetavec = cop@theta, i10 = 1:10,
                   nRnd = 50, t01 = (1:63)/64, # exact binary fractions
                   lambdaLvec = NA_real_, lambdaUvec = NA_real_)
{
    stopifnot(is(cop, "acopula"))
    cat0 <- function(...) cat(..., "\n", sep = "")
    theta0 <- cop@theta
    CT <- list()

    ## ==== (1) cop name ====

    cat0(sprintf("(1) copula family: %10s, theta0 = %g",
                 cop@name, theta0))

    ## ==== (2) generator ====

    ## ==== (2.1) psi and psiInv ====

    cat("\n(2) values of psi at i10:\n")
    CT <- c(CT, list(psi = system.time(p.i <- cop@psi(i10,theta = theta0))))
    print(p.i)
    checkifnot(identical(numeric(0), cop@psiInv(numeric(0), theta = theta0)))
    checkifnot(cop@psiInv(0, theta = theta0) == Inf)
    cat0("\nvalues of psiInv at t01:")
    CT <- c(CT, list(psiI = system.time(pi.t <-
                     cop@psiInv(t01, theta = theta0))))
    print(pi.t)
    CT[["psiI"]] <- CT[["psiI"]] +
        system.time(pi.pi <- cop@psiInv(p.i,theta = theta0))
    CT[["psi" ]] <- CT[["psi" ]] +
        system.time(p.pit <- cop@psi(pi.t, theta = theta0))
    cat0("check if psiInv(psi(i10))==i10: ", all.equal(pi.pi, i10))
    cat0("check if psi(psiInv(t01))==t01: ", all.equal(p.pit, t01))

    ## ==== (2.2) psiDAbs ====

    ## psiDAbs with degree = 1
    cat("\nvalues of psiDAbs with degree=1 at i10:\n")
    CT <- c(CT, list(psiDAbs = system.time(p.D <-
                     cop@psiDAbs(i10,theta = theta0, MC = FALSE, N = 1000))))
    print(p.D)
    cat("check if psiDAbs(Inf,theta)=0 and the class of psiDAbs(0,theta): ")
    at.0 <- cop@psiDAbs(0, theta = theta0, MC = FALSE, N = 1000)
    stopifnot(cop@psiDAbs(Inf, theta = theta0, MC = FALSE, N = 1000)==0, 
              is.numeric(at.0), !is.nan(at.0))
    cat0("[Ok]")
    ## psiDAbs with degree = 10
    cat("\nvalues of psiDAbs with degree=10 at i10:\n")
    CT <- c(CT, list(psiDAbs = system.time(p.D <- cop@psiDAbs(i10,theta = theta0,
                     degree = 10, MC = FALSE, N = 1000))))
    print(p.D)
    cat("check if psiDAbs(Inf,theta,degree=10)=0 and the class of psiDAbs(0,theta,degree=10): ")
    at.0 <- cop@psiDAbs(0, theta = theta0, degree = 10, MC = FALSE, N = 1000)
    stopifnot(cop@psiDAbs(Inf, theta = theta0, degree = 10, MC = FALSE, N = 1000) == 0,
              is.numeric(at.0), !is.nan(at.0))
    cat0("[Ok]")
    ## psiDAbs with degree = 10 and MC
    cat("\nvalues of psiDAbs with degree=10 and MC at i10:\n")
    CT <- c(CT, list(psiDAbs = system.time(p.D <- cop@psiDAbs(i10,theta = theta0,
                     degree = 10, MC = TRUE, N = 1000))))
    print(p.D)
    cat("check if all values are nonnegative\n")
    stopifnot(all(p.D >= 0))
    cat("check if psiDAbs(Inf,theta,degree=10,MC=TRUE,N=1000)=0 and the class of psiDAbs(0,theta,degree=10,MC=TRUE,N=1000): ")
    at.0 <- cop@psiDAbs(0, theta = theta0, degree = 10, MC = TRUE, N = 1000)
    stopifnot(cop@psiDAbs(Inf, theta = theta0, degree = 10, MC = TRUE, N = 1000)==0,
              is.numeric(at.0), !is.nan(at.0))
    cat0("[Ok]")

    ## ==== (3) parameter interval ====

    cat("\n(3) parameter interval:\n")
    print(cop@paraInterval)
    cat0("theta1=",theta1)
    cat0("nesting condition for theta0 and theta1 fulfilled: ",
         cop@nestConstr(theta0,theta1))

    ## ==== (4) V0, dV0, V01, dV01 ====

    ## V0
    CT <- c(CT, list(V0 = system.time(V0 <- cop@V0(nRnd,theta0))))
    cat0("\n(4) ",nRnd," generated V0's:")
    print(summary(V0))
    ## dV0
    cat("\nvalues of dV0 at i10:\n")
    CT <- c(CT, list(dV0 = system.time(dV0.i <- cop@dV0(i10,theta0))))
    print(dV0.i)
    ## V01
    CT <- c(CT, list(V01 = system.time(V01 <- cop@V01(V0,theta0,theta1))))
    cat0("\n",nRnd," generated V01's:")
    print(summary(V01))
    nt <- length(thetavec)
    ## dV01
    cat("\nvalues of dV01 at i10:\n")
    CT <- c(CT, list(dV01 = system.time(dV01.i <- cop@dV01(i10,V0=1,theta0=theta0,
                     theta1=theta1))))
    print(dV01.i)

    ## ==== (5) cCdf ====

    cat("\n(5) values of cCdf(v,0.5) for v=t01:\n")
    CT <- c(CT, list(cCdf = system.time(cCdf. <- cop@cCdf(t01,0.5,theta=theta0))))
    stopifnot(length(cCdf.) == length(t01), 0 <= cCdf.,
              TRUE ## FIXME! cCdf. <= 1
              )
    print(cCdf.)

    ## ==== (6) dAc (log = TRUE) ====

    u <- matrix(runif(1000),ncol=20)
    
    ## d = 2
    cat("\n(6) check dAc for u being a random (50x2)-matrix:\n")
    CT <- c(CT, list(dAc = system.time(lD <- cop@dAc(u[,1:2],
                     theta=theta0, MC = FALSE, N = 1000, log = TRUE))))
    stopifnot(all(is.numeric(lD) & !is.nan(lD) & is.finite(lD)))
    cat0("[Ok]")
    ## d > 2
    cat("\n check dAc for u being a random (50x20)-matrix:\n")
    CT <- c(CT, list(dAc = system.time(lD <- cop@dAc(u,theta=theta0,
                     MC = FALSE, N = 1000, log = TRUE))))
    stopifnot(all(is.numeric(lD) & !is.nan(lD) & is.finite(lD)))	
    cat0("[Ok]")
    ## d > 2, MC = TRUE
    cat("\n check dAc with MC for u being a random (50x20)-matrix:\n")
    CT <- c(CT, list(dAc = system.time(lD <- cop@dAc(u,theta=theta0, 
                     MC = TRUE, N = 1000, log = TRUE))))
    stopifnot(all(is.numeric(lD) & !is.nan(lD) & is.finite(lD)))	
    cat0("[Ok]")
    
    ## ==== (7) K ====

    check.K.t01 <- function(K)
	stopifnot(is.numeric(K), length(K) == length(t01),
		  0 <= K, K <= 1, diff(K) >= 0)
    ## K for d = 2
    cat("\n(7) values of K for d = 2 at t01:\n")
    CT <- c(CT, list(K = system.time(K. <- cop@K(t01,theta = theta0, d = 2))))
    check.K.t01( print(K.) )
    cat("check if K(0)=0 and K(1)=1: ")
    stopifnot(cop@K(0, theta = theta0, d = 2, MC = FALSE, N = 1000)==0,
              cop@K(1, theta = theta0, d = 2, MC = FALSE, N = 1000)==1)
    cat0("[Ok]")
    ## K for d = 10
    cat("\nvalues of K for d = 10 at t01:\n")
    CT <- c(CT, list(K = system.time(K. <- cop@K(t01,theta = theta0, d = 10, 
                     MC = FALSE, N = 1000))))
    check.K.t01( print(K.) )
    cat("check if K(0)=0 and K(1)=1: ")
    stopifnot(cop@K(0, theta = theta0, d = 10, MC = FALSE, N = 1000)==0,
              cop@K(1, theta = theta0, d = 10, MC = FALSE, N = 1000)==1)
    cat0("[Ok]")
    ## K for d = 10 and MC
    cat("\nvalues of K for d = 10 and MC at t01:\n")
    CT <- c(CT, list(K = system.time(K. <- cop@K(t01,theta = theta0, d = 10, 
                     MC = TRUE, N = 1000))))
    check.K.t01( print(K.) )
    cat("check if K(0)=0 and K(1)=1: ")
    stopifnot(cop@K(0, theta = theta0, d = 10, MC = TRUE, N = 1000)==0, 
              cop@K(1, theta = theta0, d = 10, MC = TRUE, N = 1000)==1)
    cat0("[Ok]")

    ## ==== (8) tau, tauInv ====

    cat("\n(8) tau at thetavec:\n")
    CT <- c(CT, list(tau = system.time(ta <- cop@tau(thetavec))))
    print(ta)
    CT <- c(CT, list(tauI = system.time(ta.I <- cop@tauInv(ta))))
    cat0("check if tauInv(tau(thetavec))==thetavec: ",
         all.equal(ta.I, thetavec))
    lambdaLvec <- rep(as.double(lambdaLvec), length.out= nt)
    lambdaUvec <- rep(as.double(lambdaUvec), length.out= nt)

    ## ==== (9) lambdaL, lambdaLInv ====

    cat("\n(9) lambdaL at thetavec:\n")
    CT <- c(CT, list(lambdaL = system.time(lT <- cop@lambdaL(thetavec))))
    CT <- c(CT, list(lT.I = system.time(lT.I <- cop@lambdaLInv(lT))))
    print(lT)
    cat0("check if lambdaLInv(lambdaL(thetavec))==lambdaLvec: ",
         all.equal(lT.I, lambdaLvec))

    ## ==== (10) lambdaU, lambdaUInv ====

    cat("\n(10) lambdaU at thetavec:\n")
    CT <- c(CT, list(lambdaU = system.time(uT <- cop@lambdaU(thetavec))))
    CT <- c(CT, list(uT.I = system.time(uT.I <- cop@lambdaUInv(uT))))
    print(uT)
    cat0("check if lambdaUInv(lambdaU(thetavec))==lambdaUvec: ",
         all.equal(uT.I, lambdaUvec))
    class(CT) <- "proc_time_list"
    CT
}

##' print() method for the tstCop() results
print.proc_time_list <- function (x, ...) {
    stopifnot(is.list(x), !is.null(nx <- names(x)))
    for(nm in nx)
        if(!all(x[[nm]] == 0, na.rm=TRUE)) {
            cat(nm,":\n"); print(x[[nm]], ...)
        }
    invisible(x)
}

## ==== copAMH ===============================================================

myAMH <- setTheta(copAMH, 0.7135001)
thetavec <- c(0.1,0.3,0.5,0.7,0.9)
tstCop(myAMH, 0.9429679, thetavec = thetavec)

## ==== copClayton ===========================================================

myClayton <- setTheta(copClayton, 0.5)
thetavec <- c(0.5,1,2,5,10)
tstCop(myClayton, 2, thetavec, lambdaL = thetavec, lambdaU = NA)

## ==== copFrank =============================================================

myFrank <- setTheta(copFrank, 1.860884)
thetavec <- c(0.5,1,2,5,10)
tstCop(myFrank, 5.736283, thetavec)

## with a slightly more extensive test:
tau.th <- c(0.055417, 0.11002, 0.21389, 0.4567, 0.66578)
tau.F <- myFrank@tau(thetavec)
stopifnot(all.equal(tau.th, tau.F, tol = 0.0001),
          all.equal(.9999, copFrank@tau(copFrank@tauInv(0.9999))),
          all.equal(myFrank@tauInv(tau.F, tol = 1e-14), thetavec, tol=1e-11))


## ==== copGumbel ============================================================

myGumbel <- setTheta(copGumbel, 1.25)
thetavec <- c(1,2,4,6,10)
tstCop(myGumbel,2, thetavec, lambdaL = NA, lambdaU = thetavec)

## ==== copJoe ===============================================================

myJoe <- setTheta(copJoe, 1.25)
thetavec <- c(1.1,2,4,6,10)
tstCop(myJoe, 2, thetavec, lambdaL = NA, lambdaU = thetavec)
