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

## A two-dimensional example: a "rotated" Clayton copula
rc3 <- rotCopula(claytonCopula(3), flip = c(TRUE, FALSE))

if(!dev.interactive(orNone=TRUE)) pdf("rotCopula-tst.pdf")
contourplot2(rc3, dCopula, nlevels = if(doExtras) 32 else 8)
contourplot2(rc3, pCopula, nlevels = if(doExtras) 32 else 8)

stopifnot(
    all.equal(rho(rc3), -0.78649216),
    all.equal(tau(rc3), -0.6))

n <- if(doExtras) 1000 else 64 # for speed
set.seed(19)
u <- rCopula(n, rc3)
(rho.n <- cor(u[,1], u[,2], method = "spearman"))
(tau.n <- cor(u[,1], u[,2], method = "kendall" ))

rC <- rotCopula(claytonCopula(), flip = c(TRUE, FALSE))
stopifnot(
    all.equal(iRho(rC, rho.n), if(doExtras) 2.877466 else 2.166027, tol=1e-6),
    all.equal(iTau(rC, tau.n), if(doExtras) 2.908705 else 2.056338, tol=1e-6))

## Fitting
system.time(f.iR  <- fitCopula(rC, pobs(u), method = "irho")); summary(f.iR)
system.time(f.iT  <- fitCopula(rC, pobs(u), method = "itau")); summary(f.iT)
system.time(f.mpl <- fitCopula(rC, pobs(u), method = "mpl" )); summary(f.mpl)

## Goodness-of-fit testing
if(doExtras >= 2)# [slow!]
    print(gC <- gofCopula(rC, u))

(gCmult <- gofCopula(rC, u, sim = "mult"))


### A four-dimensional example: a "rotated" Frank copula ---------------------
rF4 <- rotCopula(frankCopula(dim = 4), flip = c(TRUE, FALSE, TRUE, FALSE))
rF10.d4 <- setTheta(rF4, 10)

n <- if(doExtras) 1000 else 64 # for speed
set.seed(2209)
u <- rCopula(n, rF10.d4)
splom2(u)

pCopula(   c(0.6,0.7,0.6,0.8), rF10.d4) # [NaN warning - fix?]
C.n(matrix(c(0.6,0.7,0.6,0.8), 1, 4), u)

## Fitting: itau and irho should not be used (FIXME?)
f.f4 <- fitCopula(rF4, pobs(u))
summary(f.f4)

## Goodness-of-fit testing
if(doExtras)
    gofCopula(rF4, pobs(u))

gofCopula(rF4, pobs(u), sim = "mult")

## mvdc()  ex. from Johanna & Marius
rc <- rotCopula(claytonCopula(1.719))
mvcl <- mvdc(rc, c("norm", "norm"),
             list(list(mean = 0, sd = 3), list(mean = 0)))
## gave an error about @dimension previously
