## Copyright (C) 2016 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
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

## Tests for Explicit Copulas in combination with
## rotCopula, mixCopula, and khoudrajiCopula

library(copula)
isExplicit <- copula:::isExplicit
(doExtras <- copula:::doExtras())
set.seed(123)


mC <- mixCopula(list(gumbelCopula(2.5, dim = 2),
                     claytonCopula(pi, dim = 2),
                     indepCopula(dim = 2)),
                fixParam(c(2,2,4)/8, c(TRUE, TRUE, TRUE)))
mC

n <- 100
u <- as.matrix(expand.grid((1:9)/10, (1:9)/10))

## verify density and df in comparison with expression based evaluation
stopifnot(all.equal(dCopula(u, mC), copula:::dExplicitCopula.algr(u, mC)))
stopifnot(all.equal(pCopula(u, mC), copula:::pExplicitCopula.algr(u, mC)))

## derivatives
derExp <- cbind(copula:::dCdu(mC, u),    copula:::dCdtheta(mC, u),
                copula:::dlogcdu(mC, u), copula:::dlogcdtheta(mC, u))

derNum <- cbind(copula:::dCduCopulaNum(mC, u),    copula:::dCdthetaCopulaNum(mC, u),
                copula:::dlogcduCopulaNum(mC, u), copula:::dlogcdthetaCopulaNum(mC, u))

## FIXME: dlogcdu and dlogcdtheta do not match the numerical derivatives!
## Could the numerical derivates be wrong??? Why???
sapply(1:ncol(derExp), function(i) all.equal(derExp[,i], derNum[,i]))

## benchmark
## library(microbenchmark)
## microbenchmark(dCopula(u, mC), copula:::dExplicitCopula.algr(u, mC))
## microbenchmark(pCopula(u, mC), copula:::pExplicitCopula.algr(u, mC))


## rotate it
mC.surv <- rotCopula(mC)
isExplicit(mC.surv)

stopifnot(all.equal(dCopula(u, mC.surv), dCopula(1 - u, mC)))
## rotate it back
stopifnot(all.equal(dCopula(u, rotCopula(mC.surv)), dCopula(u, mC)))

## derivatives
derExp <- cbind(copula:::dCdu(mC.surv, u),    copula:::dCdtheta(mC.surv, u),
                copula:::dlogcdu(mC.surv, u), copula:::dlogcdtheta(mC.surv, u))

derNum <- cbind(copula:::dCduCopulaNum(mC.surv, u),    copula:::dCdthetaCopulaNum(mC.surv, u),
                copula:::dlogcduCopulaNum(mC.surv, u), copula:::dlogcdthetaCopulaNum(mC.surv, u))

sapply(1:ncol(derExp), function(i) all.equal(derExp[,i], derNum[,i]))

## nest the survival copula in a khoudraji Copula
k.mC.g <- khoudrajiCopula(mC.surv, gumbelCopula(3, dim = 2), c(.2, .9))
isExplicit(k.mC.g)
k.mC.g
U <- rCopula(100000, k.mC.g)

## check cdf with C.n
## cdf versus C.n
stopifnot(max(abs(pCopula(u, k.mC.g) - C.n(u, U))) < 0.002)

## pdf versus kde2d
require(MASS)
kde <- kde2d(U[,1], U[,2], n = 9, lims = c(0.1, 0.9, 0.1, 0.9))
max(abs(dCopula(u, k.mC.g) / c(kde$z) - 1)) ## relative difference < 0.185

## detivatives
derExp <- cbind(copula:::dCdu(k.mC.g, u),    copula:::dCdtheta(k.mC.g, u),
                copula:::dlogcdu(k.mC.g, u), copula:::dlogcdtheta(k.mC.g, u))

derNum <- cbind(copula:::dCduCopulaNum(k.mC.g, u),    copula:::dCdthetaCopulaNum(k.mC.g, u),
                copula:::dlogcduCopulaNum(k.mC.g, u), copula:::dlogcdthetaCopulaNum(k.mC.g, u))

sapply(1:ncol(derExp), function(i) all.equal(derExp[,i], derNum[,i]))

## nest k.mC.g and mC in a mixture copula
m.k.m <- mixCopula(list(mC, k.mC.g), c(.5, .5))
m.k.m
U <- rCopula(10000, m.k.m)
stopifnot(max(abs(pCopula(u, m.k.m) - C.n(u, U))) < 0.02)

## a monster from m.k.m and mC.surv in khoudraji
monster <- khoudrajiCopula(m.k.m, mC.surv, c(.2, .8))
monster
U <- rCopula(10000, monster)
stopifnot(max(abs(pCopula(u, monster) - C.n(u, U))) < 0.007)

## detivatives
derExp <- cbind(copula:::dCdu(monster, u),    copula:::dCdtheta(monster, u),
                copula:::dlogcdu(monster, u), copula:::dlogcdtheta(monster, u))

derNum <- cbind(copula:::dCduCopulaNum(monster, u),    copula:::dCdthetaCopulaNum(monster, u),
                copula:::dlogcduCopulaNum(monster, u), copula:::dlogcdthetaCopulaNum(monster, u))

sapply(1:ncol(derExp), function(i) all.equal(derExp[,i], derNum[,i]))

## rotate the monster
rM <- rotCopula(monster, flip=c(TRUE, FALSE))
isExplicit(rM)
rM

U <- rCopula(10000, rM)
max(abs(pCopula(u, rM) - C.n(u, U))) # < 0.005
stopifnot(identical(dCopula(u, rM), dCopula(cbind(1 - u[,1], u[,2]), monster)))
## detivatives
derExp <- cbind(copula:::dCdu(rM, u),    copula:::dCdtheta(rM, u),
                copula:::dlogcdu(rM, u), copula:::dlogcdtheta(rM, u))

derNum <- cbind(copula:::dCduCopulaNum(rM, u),    copula:::dCdthetaCopulaNum(rM, u),
                copula:::dlogcduCopulaNum(rM, u), copula:::dlogcdthetaCopulaNum(rM, u))

sapply(1:ncol(derExp), function(i) all.equal(derExp[,i], derNum[,i]))

##########################################################
## joeCopula
jC <- joeCopula(4, dim = 2)

derNum <- cbind
## pdf/cdf
stopifnot(all.equal(dCopula(u, jC), copula:::dExplicitCopula.algr(u, jC)))
stopifnot(all.equal(pCopula(u, jC), copula:::pExplicitCopula.algr(u, jC)))

## derivatives
derExp <- cbind(copula:::dCdu(jC, u),    copula:::dCdtheta(jC, u),
                copula:::dlogcdu(jC, u), copula:::dlogcdtheta(jC, u))

derNum <- cbind(copula:::dCduCopulaNum(jC, u),    copula:::dCdthetaCopulaNum(jC, u),
                copula:::dlogcduCopulaNum(jC, u), copula:::dlogcdthetaCopulaNum(jC, u))

sapply(1:ncol(derExp), function(i) all.equal(derExp[,i], derNum[,i]))

## rotate it
rJ <- rotCopula(jC, flip=c(TRUE, FALSE))

derExp <- cbind(copula:::dCdu(rJ, u),    copula:::dCdtheta(rJ, u),
                copula:::dlogcdu(rJ, u), copula:::dlogcdtheta(rJ, u))

derNum <- cbind(copula:::dCduCopulaNum(rJ, u),    copula:::dCdthetaCopulaNum(rJ, u),
                copula:::dlogcduCopulaNum(rJ, u), copula:::dlogcdthetaCopulaNum(rJ, u))

sapply(1:ncol(derExp), function(i) all.equal(derExp[,i], derNum[,i]))

## mix it with k.mc.g
hiro <- mixCopula(list(jC, k.mC.g), c(.5, .5))
hiro

##########################################################
