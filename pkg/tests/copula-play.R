library(nacopula)

myCop <- copAMH
theta(myCop) <- 0.5

cop2 <- setTheta(copAMH, value = 0.5) # is maybe more natural
stopifnot( identical(cop2, myCop) )

setGeneric("psi", function(cop) standardGeneric("psi"))
setMethod(psi, "ACopula",
          function(cop) {function(t) cop@psi(t, theta = cop@theta)})
psi(myCop) ## is a *function*
psi(myCop)(0:4)
curve(psi(myCop)(x), 0, 4)
## but this can also be done directly [ => same curve "on top" :]
curve(myCop@psi(x, theta = myCop@theta),  0, 4, col = 2, add=TRUE)

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
if(FALSE) ## FIXME
p.Tau(copJoe)# another error

#====test setup====

set.seed(1)

#====test function====

testfunction=function(mycop,theta1,thetavec,lTDCvec,uTDCvec){
  testvec=1:10
  testvec01=seq(0.1,0.9,by=0.1)
  nullvec=numeric(0)  
  n=20
  cat("(1) copula family: ",mycop@name,"\n",sep="")
  theta0=mycop@theta
  cat("theta0: ",theta0,"\n",sep="")
  cat("\n(2) values of psi at testvec:\n")
  print(mycop@psi(testvec,theta=theta0))
  cat("value of psiInv at nullvec:\n")
  print(mycop@psiInv(nullvec,theta=theta0))
  cat("values of psiInv at testvec01:\n")
  print(mycop@psiInv(testvec01,theta=theta0))
  cat("check if psiInv(psi(testvec))==testvec: ",all.equal(mycop@psiInv(mycop@psi(testvec,theta=theta0),theta=theta0),testvec),"\n",sep="")
  cat("check if psi(psiInv(testvec01))==testvec01: ",all.equal(mycop@psi(mycop@psiInv(testvec01,theta=theta0),theta=theta0),testvec01),"\n",sep="")
  cat("\n(3) parameter interval:\n")
  print(mycop@paraInterval)
  cat("nesting condition for theta0 and theta1 fulfilled: ",mycop@nestConstr(theta0,theta1),"\n",sep="")
  V0=mycop@V0(n,theta0)
  cat("\n(4) ",n," generated V0's:\n",sep="")
  print(V0)
  # V01=mycop@V01(V0,theta0,theta1)#todo: too slow for frank and joe
  # cat(n," generated V01's:\n",sep="")
  # print(V01)
  cat("\n(5) tau at thetavec:\n")
  print(mycop@tau(thetavec))
  cat("check if tauInv(tau(thetavec))==thetavec: ",all.equal(mycop@tauInv(mycop@tau(thetavec)),thetavec),"\n",sep="")
  cat("\n(6) lTDC at thetavec:\n")
  print(mycop@lTDC(thetavec))
  cat("check if lTDCInv(lTDC(thetavec))==lTDCvec: ",all.equal(mycop@lTDCInv(mycop@lTDC(thetavec)),lTDCvec),"\n",sep="")
  cat("\n(7) uTDC at thetavec:\n")
  print(mycop@uTDC(thetavec))
  cat("check if uTDCInv(uTDC(thetavec))==uTDCvec: ",all.equal(mycop@uTDCInv(mycop@uTDC(thetavec)),uTDCvec),"\n",sep="")
}

#====copAMH====

#define object
myAMH <- copAMH
theta(myAMH) <- 0.7135001
thetavec=c(0.1,0.3,0.5,0.7,0.9)

#call test function
testfunction(myAMH,0.9429679,thetavec,rep(NA*0,length(thetavec)),rep(NA*0,length(thetavec)))

#====copClayton====

#define object
myClayton <- copClayton
theta(myClayton) <- 0.5
thetavec=c(0.5,1,2,5,10)

#call test function
testfunction(myClayton,2,thetavec,thetavec,rep(NA*0,length(thetavec)))

#====copFrank===

#define object
myFrank <- copFrank
theta(myFrank) <- 1.860884
thetavec=c(0.5,1,2,5,10)

all.equal(myFrank@tauInv(myFrank@tau(thetavec)),thetavec)

#call test function
testfunction(myFrank,5.736283,thetavec,rep(NA*0,length(thetavec)),rep(NA*0,length(thetavec)))

#====copGumbel===

#define object
myGumbel <- copGumbel
theta(myGumbel) <- 1.25
thetavec=c(1,2,4,6,10)

#call test function
testfunction(myGumbel,2,thetavec,rep(NA*0,length(thetavec)),thetavec)

#====copJoe===

#define object
myJoe <- copJoe
theta(myJoe) <- 1.25
thetavec=c(1.1,2,4,6,10)

#call test function
testfunction(myJoe,2,thetavec,rep(NA*0,length(thetavec)),thetavec)

#====test value function==== 
# 
# myClayton <- copClayton
# theta(myClayton) <- 0.5
# x=c()
# value()


