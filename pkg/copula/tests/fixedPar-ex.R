## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

require(copula)

(doExtras <- copula:::doExtras())

if (doExtras)
{

### TEST FITTING ##########################################################

    n <- 100
    ## with normal copulas
    nc3  <- normalCopula(dim = 3, fixParam(c(.6,.3,.2), c(TRUE, FALSE, FALSE)),
                         dispstr = "un")
    nc3@parameters

    set.seed(4521)
    u <- pobs(rCopula(n, nc3))
    fitCopula(nc3, data = u)
    fitCopula(nc3, data = u, method = "itau")
    fitCopula(nc3, data = u, method = "irho")

    nc2  <- normalCopula(dim = 3, fixParam(c(.6,.3,.2), c(TRUE, TRUE, FALSE)),
                         dispstr = "un")
    nc2@parameters

    fitCopula(nc2, data = u)
    fitCopula(nc2, data = u, method = "itau")
    fitCopula(nc2, data = u, method = "irho")

    ## with t copulas (df.fixed = FALSE)
    tc3df  <- tCopula(dim = 3, fixParam(c(.6,.3,.2), c(TRUE, FALSE, FALSE)),
                    dispstr = "un")
    tc3df@parameters

    set.seed(4521)
    u <- pobs(rCopula(n, tc3df))
    fitCopula(tc3df, data = u)
    fitCopula(tc3df, data = u, method = "itau")
    ## fitCopula(tc3df, data = u, method = "irho")

    tc2df  <- tCopula(dim = 3, fixParam(c(.6,.3,.2), c(TRUE, TRUE, FALSE)),
                    dispstr = "un")
    tc2df@parameters

    fitCopula(tc2df, data = u)
    fitCopula(tc2df, data = u, method = "itau")
    ## fitCopula(tc2df, data = u, method = "irho")

    ## with t copulas (df.fixed = TRUE)
    tc3  <- tCopula(dim = 3, fixParam(c(.6,.3,.2), c(TRUE, FALSE, FALSE)),
                    dispstr = "un", df.fixed = TRUE)
    tc3@parameters

    fitCopula(tc3, data = u)
    fitCopula(tc3, data = u, method = "itau")
    ## fitCopula(tc3, data = u, method = "irho")

    tc2  <- tCopula(dim = 3, fixParam(c(.6,.3,.2), c(TRUE, TRUE, FALSE)),
                    dispstr = "un", df.fixed = TRUE)
    tc2@parameters

    fitCopula(tc2, data = u)
    fitCopula(tc2, data = u, method = "itau")
    ##fitCopula(tc2, data = u, method = "irho")

### TEST dC-dc functions #####################################################

    ## d*du functions should return the same result as when unfixed
    ## d*dtheta functions should return "columns" corresponding to free params
    testdCdc <- function(cop, v, cop.unfixed) {
        fixed <- attr(cop@parameters, "fixed")
        if (.hasSlot(cop, "df")) fixed <- fixed[-length(fixed)]
        stopifnot(all.equal(copula:::dCdu(cop, v), copula:::dCdu(cop.unfixed, v)),
                  all.equal(copula:::dCdtheta(cop, v),
                            copula:::dCdtheta(cop.unfixed, v)[, !fixed, drop = FALSE]),
                  all.equal(copula:::dlogcdu(cop, v), copula:::dlogcdu(cop.unfixed, v)),
                  all.equal(copula:::dlogcdtheta(cop, v),
                            copula:::dlogcdtheta(cop.unfixed, v)[, !fixed, drop = FALSE]))
    }

    ## random points in unit cube
    set.seed(7615)
    v <- matrix(runif(15), 5, 3)


    ## normal
    nc.unfixed  <- normalCopula(dim = 3, c(.6,.3,.2), dispstr = "un")
    testdCdc(nc3, v, nc.unfixed)
    testdCdc(nc2, v, nc.unfixed)

    ## t with df.fixed = TRUE
    tc.unfixed  <- tCopula(dim = 3, c(.6,.3,.2), dispstr = "un", df.fixed = TRUE)
    testdCdc(tc3, v, tc.unfixed)
    testdCdc(tc2, v, tc.unfixed)

    ## Compare true and numerical derivatives
    comparederiv <- function(cop, u) {

        c(dCdu = max(abs((copula:::dCdu(cop, u) -
                          copula:::dCduCopulaNum(cop, u)))),
          dCdtheta = max(abs(copula:::dCdtheta(cop, u) -
                             copula:::dCdthetaCopulaNum(cop, u))),
          dlogcdu = max(abs(copula:::dlogcdu(cop, u) -
                            copula:::dlogcduCopulaNum(cop, u))),
          dlogcdtheta = max(abs(copula:::dlogcdtheta(cop, u) -
                                copula:::dlogcdthetaCopulaNum(cop, u))))
    }
    comparederiv(nc3, v)
    comparederiv(nc2, v)
    comparederiv(tc3, v)
    comparederiv(tc2, v)

### Multiplier GOF #####################################################

    ## check size of mult GOF test briefly
    do1 <- function(n, cop) {
        u <- pobs(rCopula(n, cop))
        gofCopula(cop, pobs(u), sim = "mult")$p.value
    }
    M <- 10 #1000
    mean(replicate(M, do1(n, nc3)) < 0.05)
    mean(replicate(M, do1(n, nc2)) < 0.05)
    mean(replicate(M, do1(n, tc3)) < 0.05)
    mean(replicate(M, do1(n, tc2)) < 0.05)
    ## do1(n, tc3df)
    ## do1(n, tc2df)
}

