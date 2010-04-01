library(nacopula)
set.seed(1) 
n <- 2000
taueps <- 0.06

#====3d examples: generate data by hand to check if V0 and V01 are correct========================================

#correlation check function (copula sampled "by hand")
corcheck<-function(n,th0,th1,cop){
  mat <- matrix(0,nrow=n,ncol=3)
  V0time<-system.time(V0 <- cop@V0(n,th0))
  mat[,1] <- runif(n)
  V01time<-system.time(V01 <- cop@V01(V0,th0,th1))
  mat[,2] <- exp(-V0*cop@psiInv(cop@psi(rexp(n)/V01,th1),th0))
  mat[,3] <- exp(-V0*cop@psiInv(cop@psi(rexp(n)/V01,th1),th0))
  mat[,] <- cop@psi(-log(mat[,])/V0,th0)
  cormat<-cor(mat,method="kendall")
  list(V0time,V01time,cormat)
}

#create output
corcheckout<-function(corcheck,family,mat){
  cat("Run time V0 for ",family,":\n",sep="")
  print(corcheck[[1]])
  cat("\nRun time V01 for ",family,":\n",sep="")
  print(corcheck[[2]])
  cat("\nLargest distance between pairwise Kendall's taus:\n",sep="")
  print(max(abs(corcheck[[3]]-mat)))
}

#AMH
corcheckAMH<-corcheck(n,0.7135,0.9430,copAMH)
mat<-matrix(c(1,0.2,0.2,0.2,1,0.3,0.2,0.3,1),nrow=3,byrow=TRUE)#tau_{12}=tau_{13}=0.2, tau_{23}=0.3
corcheckout(corcheckAMH,"AMH",mat)
stopifnot(max(abs(corcheckAMH[[3]]-mat))<taueps)

#Clayton
corcheckClayton<-corcheck(n,0.5,2,copClayton)#tau_{12}=tau_{13}=0.2, tau_{23}=0.5
mat<-matrix(c(1,0.2,0.2,0.2,1,0.5,0.2,0.5,1),nrow=3,byrow=TRUE)#tau_{12}=tau_{13}=0.2, tau_{23}=0.5
corcheckout(corcheckClayton,"Clayton",mat)
stopifnot(max(abs(corcheckClayton[[3]]-mat))<taueps)

#Frank
corcheckFrank<-corcheck(n,1.8609,5.7363,copFrank)#tau_{12}=tau_{13}=0.2, tau_{23}=0.5
corcheckout(corcheckFrank,"Frank",mat)
stopifnot(max(abs(corcheckFrank[[3]]-mat))<taueps)

#Gumbel
corcheckGumbel<-corcheck(n,1.25,2,copGumbel)#tau_{12}=tau_{13}=0.2, tau_{23}=0.5
corcheckout(corcheckGumbel,"Gumbel",mat)
stopifnot(max(abs(corcheckGumbel[[3]]-mat))<taueps)

#Joe
corcheckJoe<-corcheck(n,1.4438,2.8562,copJoe)#tau_{12}=tau_{13}=0.2, tau_{23}=0.5
corcheckout(corcheckJoe,"Joe",mat)
stopifnot(max(abs(corcheckJoe[[3]]-mat))<taueps)

#====Examples that check value() and rn()========================================

#generate output for the examples
output=function(cormat,mat,rt){
  cat("Run time for generating ",n," vectors of variates:\n",sep="")
  print(rt)
  cat("\nLargest distance between pairwise Kendall's taus:\n",sep="")
  print(max(abs(cormat-mat)))
}

#====3d Ali-Mikhail-Haq copula example========================================

c3 <-
    new("outer_nACopula", copula = setTheta(copAMH, 0.7135),
        comp = as.integer( 1 ),
        childCops = list(new("nACopula",
                             copula = setTheta(copAMH, 0.9430),
                             comp = as.integer(c(2,3)))) # no childCops
        )

## basic check
d <- dim(c3)
stopifnot(d == 3,
	  allComp(c3) == 1:3,
	  allComp(c3@childCops[[1]]) == 2:3)

## test value()
u <- c(.3, .4, .5)
## with value function:
v <- value(c3, u)
## by hand
psi <- function(t,theta) { (1-theta)/(exp(t)-theta) }
psiInv <- function(t,theta) { log((1-theta*(1-t))/t) }
th0 <- 0.7135
th1 <- 0.9430
level1 <- psi(psiInv(u[2],th1) + psiInv(u[3],th1), th1)
level0 <- psi(psiInv(u[1],th0) + psiInv(level1, th0), th0)
stopifnot(all.equal(v, level0, tol = 1e-14))

## test rn()
rt<-system.time(rC3 <- rn(c3,n))
C3 <- cor(rC3,method="kendall")
mat<-matrix(c(1,0.2,0.2,0.2,1,0.3,0.2,0.3,1),nrow=3,byrow=TRUE)#tau_{12}=tau_{13}=0.2, tau_{23}=0.3
stopifnot(is.numeric(rC3), is.matrix(rC3),
	  dim(rC3) == c(n, 3),max(abs(C3-mat))<taueps)
output(C3,mat,rt)
# pairs(rC3, panel = function(...) {par(new=TRUE); smoothScatter(...)})

#====2d Clayton copula example========================================

c2 <-
    new("outer_nACopula", copula = setTheta(copClayton, 0.5),
        comp = as.integer(c(1,2)) # no childCops
        )
        
## basic check
d <- dim(c2)
stopifnot(d == 2,
	  allComp(c2) == 1:2)

## test value()
v <- value(c2, u = c(.3, .4))
stopifnot(all.equal(v,
 local({ u1 <- .3; u2 <- .4
        (u1^(-1/2)+u2^(-1/2)-1)^(-2)}),
                    tol = 1e-14))

## test rn()
rt<-system.time(rC2 <- rn(c2,n))
C2 <- cor(rC2,method="kendall")
mat<-matrix(c(1,0.2,0.2,1),nrow=2,byrow=TRUE)#tau_{12}=0.2
stopifnot(is.numeric(rC2), is.matrix(rC2),
	  dim(rC2) == c(n, 2), max(abs(C2-mat))<taueps)
output(C2,mat,rt)
# pairs(rC2, panel = function(...) {par(new=TRUE); smoothScatter(...)})

#====3d Clayton copula example========================================

c3 <-
    new("outer_nACopula", copula = setTheta(copClayton, 0.5),
        comp = as.integer( 1 ),
        childCops = list(new("nACopula",
                             copula = setTheta(copClayton, 2),
                             comp = as.integer(c(2,3)))) # no childCops
        )

## basic check
d <- dim(c3)
stopifnot(d == 3,
	  allComp(c3) == 1:3,
	  allComp(c3@childCops[[1]]) == 2:3)

## test value()
v <- value(c3, u = c(.3, .4, .5))
stopifnot(all.equal(v,
 local({ u1 <- .3; u2 <- .4; u3 <- .5
         1/((1/u2^2 +1/u3^2 -1)^(1/4) -1 +1/sqrt(u1))^2}),
                    tol = 1e-14))

## test rn()
rt<-system.time(rC3 <- rn(c3,n))
C3 <- cor(rC3,method="kendall")
mat<-matrix(c(1,0.2,0.2,0.2,1,0.5,0.2,0.5,1),nrow=3,byrow=TRUE)#tau_{12}=tau_{13}=0.2, tau_{23}=0.5
stopifnot(is.numeric(rC3), is.matrix(rC3),
	  dim(rC3) == c(n, 3),max(abs(C3-mat))<taueps)
output(C3,mat,rt)
# pairs(rC3, panel = function(...) {par(new=TRUE); smoothScatter(...)})

#====9d Clayton copula example========================================

c9 <- new("outer_nACopula", copula = setTheta(copClayton, 0.5),
         comp = as.integer(c(3,6,1)),
         childCops = list(new("nACopula",
                             copula = setTheta(copClayton, 2),
                             comp = as.integer(c(9,2,7,5)),
                             childCops = list(new("nACopula",
                                                  copula = setTheta(copClayton, 3),
                                                  comp = as.integer(c(8,4)) # no childCops
                                                 )
                                             )
                            )
                        )
         )

## basic check
d <- dim(c9)
stopifnot(d == 9,
    allComp(c9) == c(3,6,1,9,2,7,5,8,4),
    allComp(c9@childCops[[1]]) == c(9,2,7,5,8,4),
    allComp(c9@childCops[[1]]@childCops[[1]]) == c(8,4))

## test value()
u <- seq(0.1,0.9,by=0.1)
## with value function:
v <- value(c9, u)
## by hand
psi <- function(t,theta) { (1+t)^(-1/theta) }
psiInv <- function(t,theta) { t^(-theta) - 1 }
th0 <- 0.5
th1 <- 2
th2 <- 3
level2 <- psi(psiInv(u[8],th2) + psiInv(u[4],th2), th2)
level1 <- psi(psiInv(u[9],th1)+
              psiInv(u[2],th1)+
              psiInv(u[7],th1)+
              psiInv(u[5],th1) +
              psiInv(level2, th1), th1)
level0 <- psi(psiInv(u[3],th0)+
              psiInv(u[6],th0)+
              psiInv(u[1],th0)+
              psiInv(level1, th0), th0)
stopifnot(all.equal(v, level0, tol = 1e-14))

## test rn()
rt<-system.time(rC9 <- rn(c9,n))
C9 <- cor(rC9,method="kendall")
r1 <- c(1.0,rep(0.2,8))
r2 <- c(0.2,1.0,0.2,0.5,0.5,0.2,rep(0.5,3))
r3 <- c(0.2,0.2,1.0,rep(0.2,6))
r4 <- c(0.2,0.5,0.2,1.0,0.5,0.2,0.5,0.6,0.5)
r5 <- c(0.2,0.5,0.2,0.5,1.0,0.2,rep(0.5,3))
r6 <- c(rep(0.2,5),1.0,rep(0.2,3))
r7 <- c(0.2,0.5,0.2,0.5,0.5,0.2,1.0,0.5,0.5)
r8 <- c(0.2,0.5,0.2,0.6,0.5,0.2,0.5,1.0,0.5)
r9 <- c(0.2,0.5,0.2,0.5,0.5,0.2,0.5,0.5,1.0)
mat <- matrix(c(r1,r2,r3,r4,r5,r6,r7,r8,r9),nrow=9,byrow=TRUE)
#theoretical values: 
#C9_{11}=1.0, C9_{12}=0.2, C9_{13}=0.2, C9_{14}=0.2, C9_{15}=0.2, C9_{16}=0.2, C9_{17}=0.2, C9_{18}=0.2, C9_{19}=0.2
#             C9_{22}=1.0, C9_{23}=0.2, C9_{24}=0.5, C9_{25}=0.5, C9_{26}=0.2, C9_{27}=0.5, C9_{28}=0.5, C9_{29}=0.5
#                          C9_{33}=1.0, C9_{34}=0.2, C9_{35}=0.2, C9_{36}=0.2, C9_{37}=0.2, C9_{38}=0.2, C9_{39}=0.2
#                                       C9_{44}=1.0, C9_{45}=0.5, C9_{46}=0.2, C9_{47}=0.5, C9_{48}=0.6, C9_{49}=0.5
#                                                    C9_{55}=1.0, C9_{56}=0.2, C9_{57}=0.5, C9_{58}=0.5, C9_{59}=0.5
#                                                                 C9_{66}=1.0, C9_{67}=0.2, C9_{68}=0.2, C9_{69}=0.2
#                                                                              C9_{77}=1.0, C9_{78}=0.5, C9_{79}=0.5
#                                                                                           C9_{88}=1.0, C9_{89}=0.5
#                                                                                                        C9_{99}=1.0
stopifnot(dim(rC9) == c(n, 9),
          max(abs(C9-mat))<taueps)
output(C9,mat,rt)
# pairs(rC9, gap = .01,
# main = paste(n," vectors of a ",d,"-dimensional nested Clayton copula",sep=""))

#====Examples that check rn() according to run time========================================

#====125d Clayton copula example========================================

c125 <-
    new("outer_nACopula", copula = setTheta(copClayton, 0.5),
        comp = as.integer(),
        childCops = list(new("nACopula",
                             copula = setTheta(copClayton, 2),
                             comp = as.integer(1:10)),
                         new("nACopula",
                             copula = setTheta(copClayton, 3),
                             comp = as.integer(11:40)), 
                         new("nACopula",
                              copula = setTheta(copClayton, 2),
                              comp = as.integer(41:60)),  
                         new("nACopula",
                              copula = setTheta(copClayton, 2),
                              comp = as.integer(61:85)), 
                         new("nACopula",
                              copula = setTheta(copClayton, 3),
                              comp = as.integer(86:105)),   
                         new("nACopula",
                              copula = setTheta(copClayton, 2),
                              comp = as.integer(106:125))    
                         ) # no childCops
        )
        
## basic check
d <- dim(c125)
stopifnot(d == 125,
	  allComp(c125) == 1:125,
	  allComp(c125@childCops[[1]]) == 1:10,
	  allComp(c125@childCops[[2]]) == 11:40,
	  allComp(c125@childCops[[3]]) == 41:60,
	  allComp(c125@childCops[[4]]) == 61:85,
	  allComp(c125@childCops[[5]]) == 86:105,
	  allComp(c125@childCops[[6]]) == 106:125
	  )

## test rn()
rt<-system.time(rC125 <- rn(c125,n))
stopifnot(is.numeric(rC125), is.matrix(rC125),
	  dim(rC125) == c(n, 125))
cat("Run time for generating ",n," vectors of variates:\n",sep="")
rt
