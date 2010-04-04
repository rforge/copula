library(nacopula)
set.seed(1)
n <- 2000
N <- 50 #for the number of variates for the conducted Kolmogorov-Smirnov test
taueps <- 0.06 #for deciding if sample versions of Kendall's tau are "close enough" to population versions
doPlots <- (Sys.getenv("USER") == "maechler")

##====3d check functions=======================

##' correlation check function and run time measuring
##' @param n
##' @param th0
##' @param th1
##' @param cop
check1 <- function(n,th0,th1,cop) {
  mat <- matrix(0,nrow = n,ncol = 3)
  V0time <- system.time(V0 <- cop@V0(n,th0))
  V01time <- system.time(V01 <- cop@V01(V0,th0,th1))
  mat <- cbind(runif(n),
               exp(-V0*cop@psiInv(cop@psi(rexp(n)/V01,th1),th0)),
               exp(-V0*cop@psiInv(cop@psi(rexp(n)/V01,th1),th0)))
  mat[,] <- cop@psi(-log(mat[,])/V0,th0)
  list(V0time,V01time, cor = cor(mat,method = "kendall"))
}

##' create output
prt.tau.diff <- function(c1, c2) {
  stopifnot(is.matrix(c1))
  delta.c <- 1000 * abs(c1 - c2)
  cat(sprintf("Max & Mean distance * 1000 to true pairwise Kendall's taus: %7.1f %7.1f\n",
              max(delta.c), mean(delta.c)))
  invisible()
}
check1out <- function(x,family,trCorr) {
  cat(sprintf("Run time [ms] V0  for '%s': %7.1f\n", family, 1000*x[[1]][1]))
  cat(sprintf("Run time [ms] V01 for '%s': %7.1f\n", family, 1000*x[[2]][1]))
  prt.tau.diff(x[["cor"]], trCorr) ; cat("\n")
}

#for chisquare test: find the mass in each cube
cubemass<-function(mygridrow,cop,mesh){
  uc<-as.numeric(mygridrow)
  lc<-uc-mesh
  massincube<-value(cop,uc)-
              value(cop,c(lc[1],uc[2],uc[3]))-
              value(cop,c(uc[1],lc[2],uc[3]))-
              value(cop,c(uc[1],uc[2],lc[3]))+
              value(cop,c(lc[1],lc[2],uc[3]))+
              value(cop,c(lc[1],uc[2],lc[3]))+
              value(cop,c(uc[1],lc[2],lc[3]))-
              value(cop,lc)
  massincube
}

#for chisquare test: build a function that returns the cube number for each data point
cube<-function(data,pts){
  intervals<-as.numeric(cut(data,breaks=pts,include.lowest=TRUE,labels=1:(length(pts)-1)))
  l<-length(intervals)/3
  i1<-intervals[1:l]
  i2<-intervals[(l+1):(2*l)]
  i3<-intervals[(2*l+1):(3*l)]
  (i1-1)+(i2-1)*5+(i3-1)*5*5+1
}

#for chisquare test: define function which simulates the test statistic for 1 run
simuteststat<-function(k,n,cop,pts,cube,m,expnumofobs){
  cat("Run ",k,sep="")#user output due to possibly long run time
  data<-rn(cop,n)#generate data
  cubenumbers<-cube(data,pts)#find the cube numbers
  observationsinbin<-tabulate(cubenumbers,nbins=m)#find the number of observations in each cube
  T<-sum(((observationsinbin-expnumofobs)^2)/expnumofobs)#compute chisquare test statistic
  cat(" done\n")
  T
}

#chisquare test check function (3d only)
check2<-function(n,N,cop){

  #setup grid
  mesh<-0.2
  pts<-seq(0,1,by=mesh)#division points per dimension
  mygridall<-expand.grid(x=pts,y=pts,z=pts)#build cube upper corners
  mygrid<-mygridall[mygridall[,1]!=0 & mygridall[,2]!=0 & mygridall[,3]!=0,]#sort out those cubes with zero Lebesgue measure
  m<-length(mygrid[,1])
  row.names(mygrid)=1:m#adjust row names
  
  #determine the expected number of observations in each cube 
  masscube<-apply(mygrid,1,cubemass,cop=cop,mesh=mesh)#mass in each cube with upper corner determined by the rows of mygrid
  expnumofobs<-masscube*n#expected number of observations in each cube
  mygrid<-cbind(mygrid,expnumofobs=expnumofobs)#update grid with expected number of observations in each cube
  percentrot<-(sum(expnumofobs>=5)/m)*100#percentage of cubes that fulfill the rule of thumb
  
  #now simulate data, count observations in each cube, and compute test statistic
  T=unlist(lapply(1:N,simuteststat,n=n,cop=cop,pts=pts,cube=cube,m=m,expnumofobs=expnumofobs))
  ksresult<-ks.test(T,"pchisq",df=m-1)#compute the result of the Kolmogorov-Smirnov test based on the N realizations of the chisquare test statistics
  list(ks=ksresult,T=T,percentrot=percentrot)#return result
  
}

##====3d examples=======================

#====AMH================================

theta0<-0.7135#tau_{12}=tau_{13}=0.2, tau_{23}=0.3
theta1<-0.9430

#check 1
check1AMH <- check1(n,theta0,theta1,copAMH)
trCorr <- rbind(c(1,0.2,0.2),
                c(0.2,1,0.3),
                c(0.2,0.3,1))
check1out(check1AMH,"AMH",trCorr)
stopifnot(max(abs(check1AMH[["cor"]]-trCorr)) < taueps)

#check 2
AMH3d <-
    new("outer_nACopula", copula = setTheta(copAMH, theta0),
        comp = as.integer( 1 ),
        childCops = list(new("nACopula",
                             copula = setTheta(copAMH, theta1),
                             comp = as.integer(c(2,3)))) # no childCops
        )
resultAMH<-check2(n,N,AMH3d)
resultAMH$ks

#====Clayton================================

theta0<-0.5#tau_{12}=tau_{13}=0.2, tau_{23}=0.5
theta1<-2

#check 1
check1Clayton <- check1(n,theta0,theta1,copClayton)
trCorr <- rbind(c(1,0.2,0.2),
                c(0.2,1,0.5),
                c(0.2,0.5,1))
check1out(check1Clayton,"Clayton",trCorr)
stopifnot(max(abs(check1Clayton[["cor"]]-trCorr)) < taueps)

#check 2
Clayton3d <-
    new("outer_nACopula", copula = setTheta(copClayton, theta0),
        comp = as.integer( 1 ),
        childCops = list(new("nACopula",
                             copula = setTheta(copClayton, theta1),
                             comp = as.integer(c(2,3)))) # no childCops
        )
resultClayton<-check2(n,N,Clayton3d)
resultClayton$ks

#====Frank================================

theta0<-1.8609#tau_{12}=tau_{13}=0.2, tau_{23}=0.5
theta1<-5.7363

#check 1
check1Frank <- check1(n,theta0,theta1,copFrank)
check1out(check1Frank,"Frank",trCorr)
stopifnot(max(abs(check1Frank[["cor"]]-trCorr)) < taueps)

#check 2
Frank3d <-
    new("outer_nACopula", copula = setTheta(copFrank, theta0),
        comp = as.integer( 1 ),
        childCops = list(new("nACopula",
                             copula = setTheta(copFrank, theta1),
                             comp = as.integer(c(2,3)))) # no childCops
        )
resultFrank<-check2(n,N,Frank3d)
resultFrank$ks

#====Gumbel================================

theta0<-1.25#tau_{12}=tau_{13}=0.2, tau_{23}=0.5
theta1<-2

#check 1
check1Gumbel <- check1(n,theta0,theta1,copGumbel)
check1out(check1Gumbel,"Gumbel",trCorr)
stopifnot(max(abs(check1Gumbel[["cor"]]-trCorr)) < taueps)

#check 2
Gumbel3d <-
    new("outer_nACopula", copula = setTheta(copGumbel, theta0),
        comp = as.integer( 1 ),
        childCops = list(new("nACopula",
                             copula = setTheta(copGumbel, theta1),
                             comp = as.integer(c(2,3)))) # no childCops
        )
resultGumbel<-check2(n,N,Gumbel3d)##FIXME: log(NULL) not NULL => problem in value function!
resultGumbel$ks

#====Joe================================

theta0<-1.4438#tau_{12}=tau_{13}=0.2, tau_{23}=0.5
theta1<-2.8562

#check 1
check1Joe <- check1(n,theta0,theta1,copJoe)
check1out(check1Joe,"Joe",trCorr)
stopifnot(max(abs(check1Joe[["cor"]]-trCorr)) < taueps)

#check 2
Joe3d <-
    new("outer_nACopula", copula = setTheta(copJoe, theta0),
        comp = as.integer( 1 ),
        childCops = list(new("nACopula",
                             copula = setTheta(copJoe, theta1),
                             comp = as.integer(c(2,3)))) # no childCops
        )
resultJoe<-check2(n,N,Joe3d)
resultJoe$ks

##====Examples that check value() and rn()========================================

##generate output for the examples
prt.stats <- function(c1,c2, rt) {
  cat("Run time [ms] for generating", n,
      "vectors of variates:  ", round(1000*rt[1],1), "\n")
  prt.tau.diff(c1, c2) ; cat("\n")
}

##====3d Ali-Mikhail-Haq copula example========================================

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
rt <- system.time(rC3 <- rn(c3,n))
C3 <- cor(rC3,method = "kendall")
trCorr <- rbind(c(1,0.2,0.2),
                c(0.2,1,0.3),
                c(0.2,0.3,1))#tau_{12}=tau_{13}=0.2, tau_{23}=0.3
stopifnot(is.numeric(rC3), is.matrix(rC3),
	  dim(rC3) == c(n, 3),max(abs(C3-trCorr)) < taueps)
prt.stats(C3,trCorr,rt)
if(doPlots)
    pairs(rC3, panel = function(...) { par(new = TRUE); smoothScatter(...) })

##====2d Clayton copula example========================================

c2 <-
    new("outer_nACopula", copula = setTheta(copClayton, 0.5),
        comp = as.integer(c(1,2)) # no childCops
        )

## basic check
d <- dim(c2)
stopifnot(d ==  2,
	  allComp(c2) == 1:2)

## test value()
v <- value(c2, u = c(.3, .4))
stopifnot(all.equal(v,
 local( { u1 <- .3; u2 <- .4
        (u1^(-1/2)+u2^(-1/2)-1)^(-2) }),
                    tol = 1e-14))

## test rn()
rt <- system.time(rC2 <- rn(c2,n))
C2 <- cor(rC2,method = "kendall")
trCorr <- rbind(c(1,0.2),
                c(0.2,1))# tau_{12}=0.2
stopifnot(is.numeric(rC2), is.matrix(rC2),
	  dim(rC2) == c(n, 2), max(abs(C2-trCorr)) < taueps)
prt.stats(C2,trCorr,rt)
if(doPlots)
    smoothScatter(rC2)

##====3d Clayton copula example========================================

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
 local( { u1 <- .3; u2 <- .4; u3 <- .5
         1/((1/u2^2 +1/u3^2 -1)^(1/4) -1 +1/sqrt(u1))^2 }),
                    tol = 1e-14))

## test rn()
rt <- system.time(rC3 <- rn(c3,n))
C3 <- cor(rC3,method = "kendall")
trCorr <- matrix(c(1,0.2,0.2,0.2,1,0.5,0.2,0.5,1),nrow = 3,byrow = TRUE)#tau_{12}=tau_{13}=0.2, tau_{23}=0.5
stopifnot(is.numeric(rC3), is.matrix(rC3),
	  dim(rC3) == c(n, 3),max(abs(C3-trCorr)) < taueps)
prt.stats(C3,trCorr,rt)

if(doPlots)
    pairs(rC3, panel = function(...) { par(new = TRUE); smoothScatter(...) })

##====9d Clayton copula example========================================

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
u <- seq(0.1,0.9,by = 0.1)
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
rt <- system.time(rC9 <- rn(c9,n))
C9 <- cor(rC9,method = "kendall")
## Theoretical values:
##C_{11}= 1 , C_{12}=0.2, C_{13}=0.2, C_{14}=0.2, C_{15}=0.2, C_{16}=0.2, C_{17}=0.2, C_{18}=0.2, C_{19}=0.2
##            C_{22}= 1 , C_{23}=0.2, C_{24}=0.5, C_{25}=0.5, C_{26}=0.2, C_{27}=0.5, C_{28}=0.5, C_{29}=0.5
##                        C_{33}= 1 , C_{34}=0.2, C_{35}=0.2, C_{36}=0.2, C_{37}=0.2, C_{38}=0.2, C_{39}=0.2
##                                    C_{44}= 1 , C_{45}=0.5, C_{46}=0.2, C_{47}=0.5, C_{48}=0.6, C_{49}=0.5
##                                                C_{55}= 1 , C_{56}=0.2, C_{57}=0.5, C_{58}=0.5, C_{59}=0.5
##                                                            C_{66}= 1 , C_{67}=0.2, C_{68}=0.2, C_{69}=0.2
##                                                                        C_{77}= 1 , C_{78}=0.5, C_{79}=0.5
##                                                                                    C_{88}= 1 , C_{89}=0.5
##                                                                                                C_{99}= 1
C9.true <- rbind(c(1. ,rep(0.2,8)),
                 c(0.2,1. ,0.2,0.5,0.5,0.2, rep(0.5,3)),
                 c(0.2,0.2,1. , rep(0.2,6)),
                 c(0.2,0.5,0.2,1. ,0.5,0.2,0.5,0.6,0.5),
                 c(0.2,0.5,0.2,0.5,1. ,0.2, rep(0.5,3)),
                 c(rep(0.2,5),         1. , rep(0.2,3)),
                 c(0.2,0.5,0.2,0.5,0.5,0.2,1. ,0.5,0.5),
                 c(0.2,0.5,0.2,0.6,0.5,0.2,0.5,1. ,0.5),
                 c(0.2,0.5,0.2,0.5,0.5,0.2,0.5,0.5,1. ))
stopifnot(dim(rC9) == c(n, 9),
          max(abs(C9-C9.true)) < taueps)
prt.stats(C9,C9.true,rt)
if(doPlots && dev.interactive()) ## -> "large"
    pairs(rC9, gap = .1, pch = 20, cex = 0.2, col = rgb(.2,.1,.7, alpha = .5),
          main = paste(n," vectors of a ",d,"-dimensional nested Clayton copula",sep = ""))

##====Examples that check rn() according to run time========================================

##====125d Clayton copula example========================================

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
rt <- system.time(rC125 <- rn(c125,n))
stopifnot(is.numeric(rC125), is.matrix(rC125),
	  dim(rC125) == c(n, 125))
cat("Run time for generating ",n," vectors of variates:\n",sep = "")
rt
