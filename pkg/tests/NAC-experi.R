library(nacopula)

## Manual construction
c3 <-
    new("final_nACopula", copula = setTheta(copAMH, 1/2),
        comp = as.integer( 1 ),
        childCops = list(new("nACopula",
                             copula = setTheta(copAMH, 0.4),
                             comp = as.integer(c(2,3)))) # and no childCops
        )

d <- dim(c3)
stopifnot(d == 3,
	  allComp(c3) == 1:3,
	  allComp(c3@childCops[[1]]) == 2:3)
