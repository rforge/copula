library(nacopula)

## Manual construction
c3 <-
    new("final_nACopula", copula = setTheta(copClayton, 0.5),
        comp = as.integer( 1 ),
        childCops = list(new("nACopula",
                             copula = setTheta(copClayton, 2),
                             comp = as.integer(c(2,3)))) # and no childCops
        )

d <- dim(c3)
stopifnot(d == 3,
	  allComp(c3) == 1:3,
	  allComp(c3@childCops[[1]]) == 2:3)


##==== test value() ====
(v <- value(c3, u = c(.3, .4, .5))) # 0.1514383947453
stopifnot(all.equal(v,
 local({u1 <- .3; u2 <- .4; u3 <- .5;
        1/((1/u2^2 +1/u3^2 -1)^(1/4) -1 +1/sqrt(u1))^2}),
                    tol = 1e-14))


