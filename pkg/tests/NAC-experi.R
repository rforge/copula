library(nacopula)

## 3d Clayton copula example
c3 <-
    new("final_nACopula", copula = setTheta(copClayton, 0.5),
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
(v <- value(c3, u = c(.3, .4, .5))) # 0.1514383947453
stopifnot(all.equal(v,
 local({u1 <- .3; u2 <- .4; u3 <- .5;
        1/((1/u2^2 +1/u3^2 -1)^(1/4) -1 +1/sqrt(u1))^2}),
                    tol = 1e-14))
                         
## test rn()
rn(c3,100)              

## 9d Clayton copula example
c9 <- new("final_nACopula", copula = setTheta(copClayton, 0.5),
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
v <- value(c9, u)#with value function
#by hand
psi = function(t,theta) { (1+t)^(-1/theta) }
psiInv = function(t,theta) { t^(-theta) - 1 }
theta0=0.5
theta1=2
theta2=3
level2=psi(psiInv(u[8],theta2)+psiInv(u[4],theta2),theta2)
level1=psi(psiInv(u[9],theta1)+psiInv(u[2],theta1)+psiInv(u[7],theta1)+psiInv(u[5],theta1)+psiInv(level2,theta1),theta1)
level0=psi(psiInv(u[3],theta0)+psiInv(u[6],theta0)+psiInv(u[1],theta0)+psiInv(level1,theta0),theta0)
stopifnot(all.equal(v,level0,tol = 1e-14))                  

## test rn()
rn(c9,100)

