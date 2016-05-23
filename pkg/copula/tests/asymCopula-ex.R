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

### some constructions ###########################################################

## A Khoudraji-Clayton copula
kc <- khoudrajiCopula(copula2 = claytonCopula(6),
                      shapes = c(0.4, 0.95))
kc@parameters
contour(kc, dCopula, nlevels = 20, main = "dCopula(<khoudrajiBivCopula>)")


## True versus numerical derivatives
v <- matrix(runif(6), 3, 2)
max(abs(copula:::dCduCopulaNum(kc, v) - copula:::dCdu(kc, v)))
max(abs(copula:::dCdthetaCopulaNum(kc, v) - copula:::dCdtheta(kc, v)))

## tau, rho, lambda not supposed to work
## tau(kc)
## rho(kc)
## iTau(kc, 0.5)
## iRho(kc, 0.5)
## lambda(kc)

## A Khoudraji-Clayton copula with one fixed shape parameter
kcf <- khoudrajiCopula(copula2 = claytonCopula(6),
                       shapes = fixParam(c(0.4, 0.95), c(FALSE, TRUE)))
kcf@parameters

## True versus numerical derivatives
v <- matrix(runif(6), 3, 2)
max(abs(copula:::dCduCopulaNum(kcf, v) - copula:::dCdu(kcf, v)))
max(abs(copula:::dCdthetaCopulaNum(kcf, v) - copula:::dCdtheta(kcf, v)))

## A Khoudraji-normal-Clayton copula
knc <- khoudrajiCopula(copula1 = normalCopula(-0.7),
                       copula2 = claytonCopula(6),
                       shapes = c(0.4, 0.95))
knc@parameters
contour(knc, dCopula, nlevels = 20, main = "dCopula(<khoudrajiBivCopula>)")

## True versus numerical derivatives
max(abs(copula:::dCduCopulaNum(knc, v) - copula:::dCdu(knc, v)))
max(abs(copula:::dCdthetaCopulaNum(knc, v) - copula:::dCdtheta(knc, v)))

## A Khoudraji-normal-Clayton copula with fixed params
kncf <- khoudrajiCopula(copula1 = normalCopula(fixParam(-0.7, TRUE)),
                        copula2 = claytonCopula(6),
                        shapes = fixParam(c(0.4, 0.95), c(FALSE, TRUE)))
kncf@parameters

## True versus numerical derivatives
max(abs(copula:::dCduCopulaNum(kncf, v) - copula:::dCdu(knc, v)))
max(abs(copula:::dCdthetaCopulaNum(kncf, v) - copula:::dCdtheta(kncf, v)))


## A "nested" Khoudraji bivariate copula
kgkcf <- khoudrajiCopula(copula1 = gumbelCopula(3),
                         copula2 = kcf,
                         shapes = c(0.7, 0.25))
kgkcf@parameters
contour(kgkcf, dCopula, nlevels = 20, main = "dCopula(<khoudrajiBivCopula>)")
max(abs(copula:::dCduCopulaNum(kncf, v) - copula:::dCdu(knc, v)))
max(abs(copula:::dCdthetaCopulaNum(kncf, v) - copula:::dCdtheta(kncf, v)))


### fitting ###########################################################
n <- 300
u <- rCopula(n, kc)
plot(u)

if (doExtras)
{

    fitCopula(khoudrajiCopula(copula2 = claytonCopula()),
              start = c(1.1, 0.5, 0.5), data = pobs(u),
              optim.method = "Nelder-Mead")

    ## second shape parameter fixed to 1
    fitCopula(kcf,
              start = c(1.1, 0.5), data = pobs(u),
              optim.method = "Nelder-Mead")

    fitCopula(kcf,
              start = c(1.1, 0.5), data = pobs(u),
              optim.method = "BFGS")

    ## GOF example
    ## gofCopula(kcf, x = u, start = c(1.1, 0.5), optim.method = "BFGS")
    gofCopula(kcf, x = u, start = c(1.1, 0.5), optim.method = "BFGS", sim = "mult")

    ## check size of mult GOF test briefly
    ## do1 <- function() {
    ##     u <- rCopula(n, kc)
    ##     gofCopula(kcf, x = u, start = c(1.1, 0.5), optim.method = "BFGS",
    ##               sim = "mult")$p.value
    ## }
    ## M <- 1000
    ## res <- replicate(M, do1())
    ## mean(res < 0.05)

    ## under the alternative
    u <- rCopula(n, gumbelCopula(4))
    ## gofCopula(kcf, x = u, start = c(1.1, 0.5), optim.method = "BFGS")
    gofCopula(kcf, x = u, start = c(1.1, 0.5), optim.method = "BFGS", sim = "mult")
}

## All 'copula' subclasses
