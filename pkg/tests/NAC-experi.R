library(nacopula)

## 3d Clayton copula example
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
(v <- value(c3, u = c(.3, .4, .5))) # 0.1514384
stopifnot(all.equal(v,
 local({ u1 <- .3; u2 <- .4; u3 <- .5
         1/((1/u2^2 +1/u3^2 -1)^(1/4) -1 +1/sqrt(u1))^2}),
                    tol = 1e-14))

## test rn()
nn <- 2000
rC3 <- rn(c3,nn)
stopifnot(is.numeric(rC3), is.matrix(rC3),
	  dim(rC3) == c(nn, 3))
pairs(rC3, panel = function(...) {par(new=TRUE); smoothScatter(...)})

## 9d Clayton copula example
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
(v <- value(c9, u))# 0.0306986
## by hand
psi <- function(t,theta) { (1+t)^(-1/theta) }
psiInv <- function(t,theta) { t^(-theta) - 1 }
t0 <- 0.5
t1 <- 2
t2 <- 3
level2 <- psi(psiInv(u[8],t2) + psiInv(u[4],t2), t2)
level1 <- psi(psiInv(u[9],t1)+
              psiInv(u[2],t1)+
              psiInv(u[7],t1)+
              psiInv(u[5],t1) +
              psiInv(level2, t1), t1)
level0 <- psi(psiInv(u[3],t0)+
              psiInv(u[6],t0)+
              psiInv(u[1],t0)+
              psiInv(level1, t0), t0)
stopifnot(all.equal(v, level0, tol = 1e-14))

## test rn()
rC9 <- rn(c9,100)
stopifnot(dim(rC9) == c(100, 9))
C9 <- cor(rC9); round(C9, 3)
pairs(rC9, gap = .01,
      main = "n = 100 -- 9-dim. nested Clayton Archi.Copula")
