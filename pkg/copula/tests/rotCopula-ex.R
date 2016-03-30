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
    ## A two-dimensional example: a "rotated" Clayton copula
    rc <- rotCopula(claytonCopula(3), flip = c(TRUE, FALSE))

    ## contour(rc, dCopula, nlevels = 20)
    ## contour(rc, pCopula, nlevels = 20)
    rho(rc)
    tau(rc)

    n <- 1000
    u <- rCopula(n, rc)
    rho.n <- cor(u[,1], u[,2], method = "spearman")
    tau.n <- cor(u[,1], u[,2], method = "kendall")

    iRho(rotCopula(claytonCopula(), flip = c(TRUE, FALSE)), rho.n)
    iTau(rotCopula(claytonCopula(), flip = c(TRUE, FALSE)), tau.n)

    ## Fitting
    fitCopula(rotCopula(claytonCopula(), flip = c(TRUE, FALSE)), pobs(u),
              method = "irho")
    fitCopula(rotCopula(claytonCopula(), flip = c(TRUE, FALSE)), pobs(u),
              method = "itau")
    fitCopula(rotCopula(claytonCopula(), flip = c(TRUE, FALSE)), pobs(u),
              method = "mpl")

    ## Goodness-of-fit testing
    ## gofCopula(rotCopula(claytonCopula(), flip = c(TRUE, FALSE)), u)
    gofCopula(rotCopula(claytonCopula(), flip = c(TRUE, FALSE)), u,
              sim = "mult")

    ## A four-dimensional example: a "rotated" Frank copula
    rf <- rotCopula(frankCopula(10, dim = 4),
                    flip = c(TRUE, FALSE, TRUE, FALSE))

    n <- 1000
    u <- rCopula(n, rf)
    pairs(u)

    pCopula(c(0.6,0.7,0.6,0.8), rf)
    C.n(matrix(c(0.6,0.7,0.6,0.8), 1, 4), u)

    ## Fitting: itau and irho should not be used (FIXME?)
    fitCopula(rotCopula(frankCopula(dim=4),
                        flip = c(TRUE, FALSE, TRUE, FALSE)), pobs(u))

    ## Goodness-of-fit testing
    ## gofCopula(rotCopula(frankCopula(dim=4),
    ##                    flip = c(TRUE, FALSE, TRUE, FALSE)), pobs(u))
    gofCopula(rotCopula(frankCopula(dim=4),
                        flip = c(TRUE, FALSE, TRUE, FALSE)), pobs(u),
          sim = "mult")
}

