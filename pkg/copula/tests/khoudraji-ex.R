## Copyright (C) 2012-2017 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
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

## Khoudraji Copula

library(copula)

################################################################################
## An explicit Khoudraji copula constructed from claytonCopula and indepCoula
## with one fixed shape parameter
kcf <- khoudrajiCopula(copula2 = claytonCopula(6),
                       shapes = fixParam(c(0.4, 0.95), c(FALSE, TRUE)))
kcf
getTheta(kcf, freeOnly = FALSE, named = TRUE)
getTheta(kcf, named = TRUE)

## kcf@exprdist$cdf
stopifnot(copula:::isExplicit(kcf))

set.seed(123)
u <- rCopula(20, kcf)

## cdf/pdf from explicit expression versus from algorithmic implementation
max(abs(pCopula(u, kcf) - copula:::pKhoudrajiCopula(u, kcf)))
max(abs(dCopula(u, kcf) - copula:::dKhoudrajiBivCopula(u, kcf)))

################################################################################
## A nested Khoudraji copula from kcf and a gumbelCopula, still explicit
k_kcf_g <- khoudrajiCopula(kcf, gumbelCopula(2),
                           shapes = fixParam(c(.2, .9), c(TRUE, TRUE)))
k_kcf_g
isExplicit(k_kcf_g)
getTheta(k_kcf_g, freeOnly = FALSE, named = TRUE)
getTheta(k_kcf_g, named = TRUE)

## k_kcf_g@exprdist$cdf
set.seed(456)
U <- rCopula(100000, k_kcf_g)

## check cdf/pdf with C.n and kde2d
require(MASS)
u <- as.matrix(expand.grid((1:9)/10, (1:9)/10))
## cdf versus C.n
max(abs(pCopula(u, k_kcf_g) - C.n(u, U))) ## < 0.0001)
## pdf versus kde2d
kde <- kde2d(U[,1], U[,2], n = 9, lims = c(0.1, 0.9, 0.1, 0.9))
max(abs(dCopula(u, k_kcf_g) / c(kde$z) - 1)) ## relative difference < 0.13

################################################################################
## one more nesting: kcf and k_kcf_g
monster <- khoudrajiCopula(kcf, k_kcf_g,
                           shapes = fixParam(c(.9, .1), c(TRUE, FALSE)))
monster
isExplicit(monster)
getTheta(monster, freeOnly = FALSE, named = TRUE)
getTheta(monster, named = TRUE)

set.seed(488)
U <- rCopula(100000, monster)
## cdf versus C.n
max(abs(pCopula(u, monster) - C.n(u, U))) ## < 0.002)
## pdf versus kde2d
kde <- kde2d(U[,1], U[,2], n = 9, lims = c(0.1, 0.9, 0.1, 0.9))
max(abs(dCopula(u, monster) / c(kde$z) - 1)) ## relative difference < 0.09

################################################################################
## khoudrajiExplicitCopula with dim 3
kcd3 <- khoudrajiCopula(copula1 = indepCopula(dim=3),
                        copula2 = claytonCopula(6, dim=3),
                        shapes = c(0.4, 0.95, 0.95))
kcd3
isExplicit(kcd3)
class(kcd3) ## "khoudrajiExplicitCopula"

set.seed(1712)
U <- rCopula(100000, kcd3)
u <- matrix(runif(15), 5, 3)
## cdf versus C.n
max(abs(pCopula(u, kcd3) - C.n(u, U))) # < 0.00007
## don't know how to check pdf
(f.v <- dCopula(u, kcd3))
stopifnot(
    all.equal(f.v,
              c(1.1116773, 1.6661734, 0.1719080, 0.4981574, 9.8964415),
              tol = 1e-7)
)

kcg <- khoudrajiCopula(claytonCopula(2, dim=3), gumbelCopula(3, dim=3), c(0.2,  0.2, .8))
kcg
u <- rCopula(10, kcg)
dCopula(u, kcg)


kkcgg <- khoudrajiCopula(kcg, gumbelCopula(2, dim = 3), c(.3, .3, .9))
kkcgg
u <- rCopula(10, kkcgg)
dCopula(u, kkcgg)


kcgcd3 <- khoudrajiCopula(kcd3, kkcgg, c(.4, .4, .5))
kcgcd3
u <- rCopula(10, kcgcd3)
dCopula(u, kcgcd3)

## dCdu 
all.equal(copula:::dCdu(kcgcd3, u), copula:::dCduCopulaNum(kcgcd3, u))
## dCdtheta
all.equal(copula:::dCdtheta(kcgcd3, u), copula:::dCdthetaCopulaNum(kcgcd3, u))
## dlogcdu
all.equal(copula:::dlogcdu(kcgcd3, u), copula:::dlogcduCopulaNum(kcgcd3, u))
## dlogcdtheta
all.equal(copula:::dlogcdtheta(kcgcd3, u), copula:::dlogcdthetaCopulaNum(kcgcd3, u))

#########################################################################################
## fitting checking

do1 <- function(n, cop) {
    U <- rCopula(n, cop)
    fit <- try(fitCopula(cop, U))
    p <- copula:::nParam(cop, freeOnly=TRUE)
    if (inherits(fit, "try-error")) rep(NA_real_, 2 * p)
    else c(coef(fit), sqrt(diag(vcov(fit))))
}

sumsim <- function(sim) {
    p <- nrow(sim) / 2
    mm <- apply(sim, 1, mean, na.rm = TRUE)
    ss <- apply(sim, 1, sd, na.rm = TRUE)
    data.frame(estimate = mm[1:p], ese = ss[1:p], ase = mm[p + 1:p])
}

set.seed(123)
nrep <- 50
n <- 1000

sim <- replicate(nrep, do1(n, kcg))
sumsim(sim)
cbind(true = getTheta(kcg), sumsim(sim)) ## the average se is greater than empirical; need more digging
