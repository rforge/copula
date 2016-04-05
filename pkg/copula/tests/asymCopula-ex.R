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

if(!dev.interactive(orNone=TRUE)) pdf("asymCopula-ex.pdf")

## An asymetric Clayton copula object
ac <- asymBivCopula(copula1 = indepCopula(),
                    copula2 = claytonCopula(6),
                    shapes = c(0.4, 0.95))
contour(ac, dCopula, nlevels = 20, main = "dCopula(<asymBivCopula>)")

## true versus numerical derivatives
v <- matrix(runif(6), 3, 2)
max(abs(copula:::dCduCopulaNum(ac, v) - copula:::dCdu(ac, v)))
max(abs(copula:::dCdthetaCopulaNum(ac, v) - copula:::dCdtheta(ac, v)))

## tau, rho, lambda not supposed to work
## tau(ac)
## rho(ac)
## iTau(ac, 0.5)
## iRho(ac, 0.5)
## lambda(ac)

## fitting example
n <- 300
u <- rCopula(n, ac)
##plot(u)

if (doExtras)
{

    fitCopula(asymBivCopula(copula2 = claytonCopula()),
              start = c(1.1, 0.5, 0.5), data = pobs(u),
              optim.method="Nelder-Mead")

    ## GOF example

    ## gofCopula(asymBivCopula(copula2 = claytonCopula()), pobs(u),
    ##           start = c(1.1, 0.5, 0.5), optim.method="Nelder-Mead")

    gofCopula(asymBivCopula(copula2 = claytonCopula()), pobs(u),
              start = c(1.1, 0.5, 0.5), optim.method="Nelder-Mead", sim = "mult")

    ## check size of GOF test briefly
    ## do1 <- function() {
    ##     u <- rCopula(n, ac)
    ##     gofCopula(asymBivCopula(copula2 = claytonCopula()), pobs(u),
    ##               start = c(1.1, 0.5, 0.5), optim.method="Nelder-Mead",
    ##               sim = "mult")$p.value
    ## }
    ## M <- 100
    ## res <- replicate(M, do1())
    ## mean(res < 0.05)

    ## under the alternative
    u <- rCopula(n, gumbelCopula(4))

    ## gofCopula(asymBivCopula(copula2 = claytonCopula()), pobs(u),
    ##          start = c(1.1, 0.5, 0.5), optim.method="Nelder-Mead")

    gofCopula(asymBivCopula(copula2 = claytonCopula()), pobs(u),
              start = c(1.1, 0.5, 0.5), optim.method="Nelder-Mead", sim = "mult")

    ## a "nested" asymetric bivariate copula
    agac <- asymBivCopula(copula1 = gumbelCopula(3),
                          copula2 = ac,
                          shapes = c(0.4, 0.95))


}

## All 'copula' subclasses
