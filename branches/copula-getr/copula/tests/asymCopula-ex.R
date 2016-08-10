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
kc
contour(kc, dCopula, nlevels = 20, main = "dCopula(<khoudrajiBivCopula>)")

## check density: special case where we know the answer
kd2a <- khoudrajiCopula(copula1 = indepCopula(),
                        copula2 = claytonCopula(6),
                        shapes = c(1, 1))

kd2b <- khoudrajiCopula(copula1 = gumbelCopula(4),
                        copula2 = indepCopula(),
                        shapes = c(0, 0))

v <- matrix(runif(10), 5, 2)

stopifnot(all.equal(dCopula(v, kd2a), dCopula(v, claytonCopula(6))))
stopifnot(all.equal(dCopula(v, kd2b), dCopula(v, gumbelCopula(4))))


## True versus numerical derivatives
v <- matrix(runif(6), 3, 2)
max(abs(copula:::dCduCopulaNum(kc, v) - copula:::dCdu(kc, v)))
max(abs(copula:::dCdthetaCopulaNum(kc, v) - copula:::dCdtheta(kc, v)))

## tau, rho, lambda not supposed to work
try(tau(kc))
try(rho(kc))
try(iTau(kc, 0.5))
try(iRho(kc, 0.5))
try(lambda(kc))

## A Khoudraji-Clayton copula with one fixed shape parameter
kcf <- khoudrajiCopula(copula2 = claytonCopula(6),
                       shapes = fixParam(c(0.4, 0.95), c(FALSE, TRUE)))
kcf


## True versus numerical derivatives
v <- matrix(runif(6), 3, 2)
max(abs(copula:::dCduCopulaNum(kcf, v) - copula:::dCdu(kcf, v)))
max(abs(copula:::dCdthetaCopulaNum(kcf, v) - copula:::dCdtheta(kcf, v)))

## A Khoudraji-normal-Clayton copula
knc <- khoudrajiCopula(copula1 = normalCopula(-0.7),
                       copula2 = claytonCopula(6),
                       shapes = c(0.4, 0.95))
knc
contour(knc, dCopula, nlevels = 20, main = "dCopula(<khoudrajiBivCopula>)")

## True versus numerical derivatives
max(abs(copula:::dCduCopulaNum(knc, v) - copula:::dCdu(knc, v)))
max(abs(copula:::dCdthetaCopulaNum(knc, v) - copula:::dCdtheta(knc, v)))

## A Khoudraji-normal-Clayton copula with fixed params
kncf <- khoudrajiCopula(copula1 = normalCopula(fixParam(-0.7, TRUE)),
                        copula2 = claytonCopula(6),
                        shapes = fixParam(c(0.4, 0.95), c(FALSE, TRUE)))
kncf

## True versus numerical derivatives
max(abs(copula:::dCduCopulaNum(kncf, v) - copula:::dCdu(knc, v)))
max(abs(copula:::dCdthetaCopulaNum(kncf, v) - copula:::dCdtheta(kncf, v)))


## A "nested" Khoudraji bivariate copula
kgkcf <- khoudrajiCopula(copula1 = gumbelCopula(3),
                         copula2 = kcf,
                         shapes = c(0.7, 0.25))
kgkcf
contour(kgkcf, dCopula, nlevels = 20, main = "dCopula(<khoudrajiBivCopula>)")
max(abs(copula:::dCduCopulaNum(kgkcf, v) - copula:::dCdu(kgkcf, v)))
max(abs(copula:::dCdthetaCopulaNum(kgkcf, v) - copula:::dCdtheta(kgkcf, v)))

## A three dimensional Khoudraji-Clayton copula
kcd3 <- khoudrajiCopula(copula1 = indepCopula(dim=3),
                        copula2 = claytonCopula(6, dim=3),
                        shapes = c(0.4, 0.95, 0.95))
kcd3
class(kcd3) ## "khoudrajiExplicitCopula"
set.seed(1712)
n <- 300
u <- rCopula(n, kcd3)
v <- matrix(runif(15), 5, 3)
splom2(u)
dCopula(v, kcd3)

## check density: special case where we know the answer
kd3a <- khoudrajiCopula(copula1 = indepCopula(dim=3),
                        copula2 = claytonCopula(6, dim=3),
                        shapes = c(1, 1, 1))

kd3b <- khoudrajiCopula(copula1 = gumbelCopula(4, dim=3),
                        copula2 = indepCopula(dim=3),
                        shapes = c(0, 0, 0))

stopifnot(all.equal(dCopula(v, kd3a), dCopula(v, claytonCopula(6, dim=3))))
stopifnot(all.equal(dCopula(v, kd3b), dCopula(v, gumbelCopula(4, dim=3))))


## A four dimensional Khoudraji-Normal copula
knd4 <- khoudrajiCopula(copula1 = indepCopula(dim=4),
                        copula2 = normalCopula(.9, dim=4),
                        shapes = c(0.4, 0.95, 0.95, 0.95))
knd4
class(knd4) ## "khoudrajiCopula"
u <- rCopula(n, knd4)
splom2(u)
v <- matrix(runif(20), 5, 4)
try(dCopula(v, knd4)) ## should fail

## a nested Khoudraji copula whose construction should result
## in a khoudrajiCopula object, but not a khoudrajiExplicitCopula object

kd3 <- khoudrajiCopula(copula1 = khoudrajiCopula(copula1 = gumbelCopula(4, dim=3),
                                                 copula2 = claytonCopula(6, dim=3),
                                                 shapes = c(.4, 0.95, 0.95)),
                       copula2 = frankCopula(10, dim=3),
                       shapes = c(.4, 0.95, 0.95))

kd3
class(kd3) # this should be a khoudrajiCopula, not a khoudrajiExplicitCopula
## as second argument copula is not symmetric (and in practice is not archmCopula)

## GETR
## A Khoudraji copula constructed from the survival Clayton
krc <- khoudrajiCopula(copula2 = rotCopula(claytonCopula(6)),
                      shapes = c(0.4, 0.95))
krc
#contour(krc, dCopula, nlevels = 20, main = "dCopula(<khoudrajiBivCopula>)")
u <- rCopula(n, krc)
splom2(u)

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
    summary(fk1)

    ## kcf : second shape parameter fixed to 0.95
    fkN <- fitCopula(kcf,
                     start = c(1.1, 0.5), data = pobs(u),
                     optim.method = "Nelder-Mead", optim.control = list(trace = TRUE))
    fkN
    summary(fkN)
    fkB <- fitCopula(kcf,
                     start = c(1.1, 0.5), data = pobs(u),
                     optim.method = "BFGS", optim.control = list(trace = TRUE))
    fkB
    summary(fkB)
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

}

if (FALSE)
{
    ##############################################################################
    ## Testing higher-dimensional density / fitting

    n <- 500

    set.seed(1251)

    ## A three dimensional Khoudraji-Clayton copula
    kcd3 <- khoudrajiCopula(copula1 = indepCopula(dim=3),
                            copula2 = claytonCopula(6, dim=3),
                            shapes = c(0.4, 0.95, 0.95))

    ## one fitting
    u <- rCopula(n, kcd3)
    fitCopula(khoudrajiCopula(copula1 = indepCopula(dim=3),
                              copula2 = claytonCopula(dim=3)),
              start = c(1.1, 0.5, 0.5, 0.5), data = pobs(u),
              optim.method = "Nelder-Mead")


    ## bias and stderr
    do1 <- function() {
        u <- rCopula(n, kcd3)
        fit <- fitCopula(khoudrajiCopula(copula1 = indepCopula(dim=3),
                                         copula2 = claytonCopula(dim=3)),
                         start = c(1.1, 0.5, 0.5, 0.5), data = pobs(u),
                         optim.method = "Nelder-Mead", estimate.variance = TRUE)
        c(fit@estimate, sqrt(diag(fit@var.est)))
    }
    res <- replicate(100, do1())
    apply(res,1,mean)
    apply(res,1,sd)

    ## A three dimensional Khoudraji-Gumbel copula
    kgd3 <- khoudrajiCopula(copula1 = indepCopula(dim=3),
                            copula2 = gumbelCopula(6, dim=3),
                            shapes = c(0.4, 0.95, 0.95))

    ## one fitting
    u <- rCopula(n, kgd3)
    fitCopula(khoudrajiCopula(copula1 = indepCopula(dim=3),
                              copula2 = gumbelCopula(dim=3)),
              start = c(1.1, 0.5, 0.5, 0.5), data = pobs(u),
              optim.method = "Nelder-Mead")


    ## bias and stderr
    do1 <- function() {
        u <- rCopula(n, kgd3)
        fit <- fitCopula(khoudrajiCopula(copula1 = indepCopula(dim=3),
                                         copula2 = gumbelCopula(dim=3)),
                         start = c(1.1, 0.5, 0.5, 0.5), data = pobs(u),
                         optim.method = "Nelder-Mead", estimate.variance = TRUE)
        c(fit@estimate, sqrt(diag(fit@var.est)))
    }
    res <- replicate(100, do1())
    apply(res,1,mean)
    apply(res,1,sd)

    ## A three dimensional Khoudraji-Frank copula
    kfd3 <- khoudrajiCopula(copula1 = indepCopula(dim=3),
                            copula2 = frankCopula(10, dim=3),
                            shapes = c(0.4, 0.95, 0.95))

    ## one fitting
    u <- rCopula(n, kfd3)
    fitCopula(khoudrajiCopula(copula1 = indepCopula(dim=3),
                              copula2 = frankCopula(dim=3)),
              start = c(1.1, 0.5, 0.5, 0.5), data = pobs(u),
              optim.method = "Nelder-Mead")


    ## bias and stderr
    do1 <- function() {
        u <- rCopula(n, kfd3)
        fit <- fitCopula(khoudrajiCopula(copula1 = indepCopula(dim=3),
                                         copula2 = frankCopula(dim=3)),
                         start = c(1.1, 0.5, 0.5, 0.5), data = pobs(u),
                         optim.method = "Nelder-Mead", estimate.variance = TRUE)
        c(fit@estimate, sqrt(diag(fit@var.est)))
    }
    res <- replicate(100, do1())
    apply(res,1,mean)
    apply(res,1,sd)
}

## All 'copula' subclasses
