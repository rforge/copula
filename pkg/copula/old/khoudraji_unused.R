## Moved from ../R/asymCopula.R -- after the better  dKhoudrajiExplicitCopula.algr  were introduced:

#### The following is no longer used #################################################################
#### They are much slower than the one above
## > microbenchmark(dKhoudrajiExplicitCopula.algr(u, kcd3), dCopula(u, kcd3))
## Unit: milliseconds
##                                    expr       min        lq      mean    median
##  dKhoudrajiExplicitCopula.algr(u, kcd3)  2.107368  2.268947  2.660072  2.471952
##                        dCopula(u, kcd3) 60.112831 64.307831 67.056729 66.658310
##         uq      max neval
##   2.569195 19.22082   100
##  68.780618 89.72618   100


##' @title Get the power set of a sequence 1:d
##' @param d the dimension
##' @return a logical matrix of 2^d rows and d columns
##' where each row indicates one subset (TRUE indicating a component is selected)
getPowerSet <- function(d) { ## TODO: impossible for large d -- need "nextSet()" there!
    as.matrix(unname(expand.grid(replicate(d, list(c(TRUE,FALSE))),
                                 KEEP.OUT.ATTRS = FALSE)))
}

##' @title The derivatives of CDF of a copula needed in Khoudraji's copula density
##' @param idx logical vector of d dimension, TRUE means derivative requested;
##' it is designed to be one row from the returned matrix from getPowerSet.
##' @param u matrix of observations with d columns at which the derivatives are needed
##' @param dg matrix of derivatives of g at u
##' @param copula a copula object
##' @param derExprs expressions of length d + 1, the first one being the cdf of the copula,
##' the second one dC / du1, the third one dC2 / (du1 du2), and so on.
##' @return a vector of the derivatives requested via idx
##' @author Jun Yan
##' __NOT_EXPORTED__ ==> not checking arguments
densDers <- function(idx, u, dg, copula, derExprs) {
    ## assuming exchangeable copula1 and copula2
    ## IK: assuming one-parameter copulas also
    dorder <- sum(idx)
    alpha <- copula@parameters[1] # possibly needed in 'derExprs' below
    d <- copula@dimension
    newidx <- c((1:d)[idx], (1:d)[!idx])
    u <- u[, newidx]
    for (i in 1:d) assign(paste0("u", i), u[,i])
    dgu <- if (sum(idx) == 0) 1 else apply(dg[,idx,drop=FALSE], 1, prod)
    c(eval(derExprs[dorder + 1])) * dgu
}


dKhoudrajiExplicitCopula <- function(u, copula, log = FALSE, ...) {
    d <- copula@dimension
    if (!is.matrix(u)) u <- matrix(u, ncol = d)
    comps <- getKhoudrajiCopulaComps(copula)
    copula1 <- comps$copula1
    copula2 <- comps$copula2
    a <- matrix(comps$shapes, nrow(u), ncol(u), byrow=TRUE)
    u1 <- u ^ (1 - a)
    u2 <- u ^ a
    dg1 <- (1 - a) * u^(-a)
    dg2 <- a * u^(a - 1)
    powerSet <- getPowerSet(d)
    dens <- 0
    for (i in 1:nrow(powerSet)) {
        idx1 <- c(powerSet[i,])
        idx2 <- c(!powerSet[i,])
        ## WARNING: Again, this works for exchangeable copula components only
        part1 <- densDers(idx1, u1, dg1, copula1, copula@derExprs1)
        part2 <- densDers(idx2, u2, dg2, copula2, copula@derExprs2)
        dens <- dens + part1 * part2
        ## print(part1); print(part2)
    }
    if(log) log(dens) else dens
}
