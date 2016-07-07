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
max(abs(copula:::dCduCopulaNum(kgkcf, v) - copula:::dCdu(kgkcf, v)))
max(abs(copula:::dCdthetaCopulaNum(kgkcf, v) - copula:::dCdtheta(kgkcf, v)))

## A three dimensional Khoudraji-Clayton copula
kcd3 <- khoudrajiCopula(copula1 = indepCopula(dim=3),
                        copula2 = claytonCopula(6, dim=3),
                        shapes = c(0.4, 0.95, 0.95))
kcd3@parameters
class(kcd3) ## "khoudrajiExplicitCopula"
set.seed(1712)
n <- 1000
u <- rCopula(n, kcd3)
v <- matrix(runif(15), 5, 3)
splom2(u)
try(dCopula(v, kcd3)) ## bugged

## A four dimensional Khoudraji-Normal copula
knd4 <- khoudrajiCopula(copula1 = indepCopula(dim=4),
                        copula2 = normalCopula(.9, dim=4),
                        shapes = c(0.4, 0.95, 0.95, 0.95))
knd4
class(knd4) ## "khoudrajiCopula"
u <- rCopula(n, knd4)
splom2(u)
try(dCopula(v, knd4)) ## not implemented


## comparing copula density between dKhoudrajiCopula and dKhoudrajiExplicitCopula
kgc    <-                 khoudrajiCopula(gumbelCopula(2), claytonCopula(6), c(.4, .95))
kgc.ex <- copula:::khoudrajiExplictCopula(gumbelCopula(2), claytonCopula(6), c(.4, .95))

u <- rCopula(20, kgc)
copula:::dKhoudrajiBivCopula(u, kgc)
copula:::dKhoudrajiExplicitCopula(u, kgc.ex)


### fitting ###########################################################
n <- 300
set.seed(17)
u <- rCopula(n, kc)
plot(u)

if (doExtras)
{
    fk1 <- fitCopula(khoudrajiCopula(copula2 = claytonCopula()),
                     start = c(1.1, 0.5, 0.5), data = pobs(u),
                     optim.method = "Nelder-Mead", optim.control = list(trace = TRUE))
    fk1

    ## kcf : second shape parameter fixed to 0.95
    fkN <- fitCopula(kcf,
                     start = c(1.1, 0.5), data = pobs(u),
                     optim.method = "Nelder-Mead", optim.control = list(trace = TRUE))
    fkN
    fkB <- fitCopula(kcf,
                     start = c(1.1, 0.5), data = pobs(u),
                     optim.method = "BFGS", optim.control = list(trace = TRUE))
    fkB
    stopifnot(
        all.equal(coef(fk1), c(c2.param = 5.42332, shape1 = 0.364467, shape2 = 0.868297),
                  tol = 1e-4), # seen 2.7e-7
        all.equal(coef(fkN), c(c2.param = 4.40525, shape1 = 0.389693), tol = 1e-4), # seen 3e-7
        all.equal(coef(fkN), coef(fkB), tol = 1e-3) # seen 1.19e-4
    )

    ## GOF example
    ## takes too long:
    ## gofCopula(kcf, x = u, start = c(1.1, 0.5), optim.method = "BFGS")
    set.seed(12)
    g.kcf <- gofCopula(kcf, x = u, start = c(1.1, 0.5), optim.method = "BFGS", sim = "mult")
    g.kcf
    stopifnot(inherits(g.kcf, "htest"),
              all.equal(g.kcf$p.value, 0.05544456, tol = 1e-3) # seen 8.2e-8
              )

    ## check size of mult GOF test briefly
    ## do1 <- function() {
    ##     u <- rCopula(n, kc)
    ##     gofCopula(kcf, x = u, start = c(1.1, 0.5), optim.method = "Nelder-Mead",
    ##               sim = "mult")$p.value
    ## }
    ## M <- 1000
    ## res <- replicate(M, do1())
    ## mean(res < 0.05)

    ## under the alternative
    set.seed(11); u <- rCopula(n, gumbelCopula(4))
    ## gofCopula(kcf, x = u, start = c(1.1, 0.5), optim.method = "BFGS")
    gA <- try(gofCopula(kcf, x = u, start = c(1.1, 0.5), optim.method = "BFGS", sim = "mult"))
    ## optim(*, "BFGS") sometimes not converging, but the above ex. does and  "Nelder-Mead" does:
    set.seed(17)
    g2 <- gofCopula(kcf, x = u, start = c(1.1, 0.5), optim.method = "Nelder-M", sim = "mult")
    stopifnot(inherits(g2, "htest"),
	      all.equal(g2$p.value, 0.0004995005, tol = 1e-4))# seen 1e-9

    ##############################################################################
    ## for JY: testing higher-dimensional density
    ## copula:::dKhoudrajiExplicitCopula(v, kcd3)
    ## copula:::dCopula(matrix(v, kcd3)) ## should the same

    ## n <- 1000
    ## do1 <- function() {
    ##     u <- rCopula(n, kcd3)
    ##     fitCopula(khoudrajiCopula(copula1 = indepCopula(dim=3),
    ##                               copula2 = claytonCopula(dim=3)),
    ##               start = c(1.1, 0.5, 0.5, 0.5), data = pobs(u),
    ##               optim.method = "Nelder-Mead")$estimate
    ## }
    ## M <- 10 ## 1000
    ## res <- replicate(M, do1())
    ## mean(res)
    ## sd(res)

}

## All 'copula' subclasses
