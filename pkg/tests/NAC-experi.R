library(nacopula)

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
(v <- value(c2, u = c(.3, .4)))
stopifnot(all.equal(v,
 local({ u1 <- .3; u2 <- .4
        (u1^(-1/2)+u2^(-1/2)-1)^(-2)}),
                    tol = 1e-14))

## test rn()
nn <- 2000
rC2 <- rn(c2,nn)
stopifnot(is.numeric(rC2), is.matrix(rC2),
	  dim(rC2) == c(nn, 2))
C2 <- cor(rC2,method="kendall")
formatC(round(C2, 3),format="f",digits=3)#theoretical value: 0.2
pairs(rC2, panel = function(...) {par(new=TRUE); smoothScatter(...)})

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
(v <- value(c3, u))
## by hand
psi <- function(t,theta) { (1-theta)/(exp(t)-theta) }
psiInv <- function(t,theta) { log((1-theta*(1-t))/t) }
th0 <- 0.7135
th1 <- 0.9430
level1 <- psi(psiInv(u[2],th1) + psiInv(u[3],th1), th1)
level0 <- psi(psiInv(u[1],th0) + psiInv(level1, th0), th0)
stopifnot(all.equal(v, level0, tol = 1e-14))

## test rn()
nn <- 2000
rC3 <- rn(c3,nn)
stopifnot(is.numeric(rC3), is.matrix(rC3),
	  dim(rC3) == c(nn, 3))
C3 <- cor(rC3,method="kendall")
formatC(round(C3, 3),format="f",digits=3)#theoretical values: C3_{12}=0.2, C3_{13}=0.2, C3_{23}=0.3
pairs(rC3, panel = function(...) {par(new=TRUE); smoothScatter(...)})

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
(v <- value(c3, u = c(.3, .4, .5)))
stopifnot(all.equal(v,
 local({ u1 <- .3; u2 <- .4; u3 <- .5
         1/((1/u2^2 +1/u3^2 -1)^(1/4) -1 +1/sqrt(u1))^2}),
                    tol = 1e-14))

## test rn()
nn <- 2000
rC3 <- rn(c3,nn)
stopifnot(is.numeric(rC3), is.matrix(rC3),
	  dim(rC3) == c(nn, 3))
C3 <- cor(rC3,method="kendall")
formatC(round(C3, 3),format="f",digits=3)#theoretical values: C3_{12}=0.2, C3_{13}=0.2, C3_{23}=0.5
pairs(rC3, panel = function(...) {par(new=TRUE); smoothScatter(...)})

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
(v <- value(c9, u))
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
rC9 <- rn(c9,100)
stopifnot(dim(rC9) == c(100, 9))
C9 <- cor(rC9,method="kendall")
formatC(round(C9, 3),format="f",digits=3)
#theoretical values: 
#C9_{12}=0.2, C9_{13}=0.2, C9_{14}=0.2, C9_{15}=0.2, C9_{16}=0.2, C9_{17}=0.2, C9_{18}=0.2, C9_{19}=0.2
#             C9_{23}=0.2, C9_{24}=0.5, C9_{25}=0.5, C9_{26}=0.2, C9_{27}=0.5, C9_{28}=0.5, C9_{29}=0.5
#                          C9_{34}=0.2, C9_{35}=0.2, C9_{36}=0.2, C9_{37}=0.2, C9_{38}=0.2, C9_{39}=0.2
#                                       C9_{45}=0.5, C9_{46}=0.2, C9_{47}=0.5, C9_{48}=0.6, C9_{49}=0.5
#                                                    C9_{56}=0.2, C9_{57}=0.5, C9_{58}=0.5, C9_{59}=0.5
#                                                                 C9_{67}=0.2, C9_{68}=0.2, C9_{69}=0.2
#                                                                              C9_{78}=0.5, C9_{79}=0.5
#                                                                                           C9_{89}=0.5
pairs(rC9, gap = .01,
      main = "n = 100 -- 9-dim. nested Clayton Archi.Copula")
