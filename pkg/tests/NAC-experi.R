library(nacopula)

## Manual construction
## --- Currently *FAILS* (in validity checking),
##     as  the subcopula has  dim(.) = 2, but contains components outside (1:2)
c3 <-
    new("nACopula", copula = setTheta(copAMH, 1/2),
        comp = as.integer( 1 ),
        childCops = list(new("nACopula",
                             copula = setTheta(copAMH, 0.4),
                             comp = as.integer(c(2,3)))) # and no childCops
        )

