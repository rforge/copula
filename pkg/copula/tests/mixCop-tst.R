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

### Finite Mixtures of Copulas:  mixCopula()
### ========================================

library(copula)
isExplicit <- copula:::isExplicit
(doExtras <- copula:::doExtras())

## if(!dev.interactive(orNone=TRUE)) pdf("mixCop-tst.pdf")

mC <- mixCopula(list(gumbelCopula(2.5, dim=3),
                     claytonCopula(pi, dim=3),
                     tCopula(0.7, dim=3)),
                c(2,2,4)/8)
mC
stopifnot(dim(mC) == 3, inherits(mC, "mixCopula"))

## mix a Gumbel with a rotated Gumbel (with equal weights 1/2):
mGG <- mixCopula(list(gumbelCopula(2), rotCopula(gumbelCopula(1.5))))
stopifnot(dim(mGG) == 2, inherits(mGG, "mixCopula"),
    all.equal(   rho(mGG), 0.57886340158, tol=1e-10) ,
    all.equal(lambda(mGG), c(lower = 0.206299474016,
                             upper = 0.292893218813), tol = 1e-10)
)
mGG # 4 parameters

set.seed(17)
uM  <- rCopula( 600, mC)
uGG <- rCopula(1000, mGG)

stopifnot(
    all.equal(
        loglikCopula(c(2.5, 1., w = c(2,4)/6), u=uGG, copula = mGG),
        37.04012, tol = 7e-7),
    all.equal(
        loglikCopula(c(2.5, pi, rho.1=0.7, df = 4, w = c(2,2,4)/8),
                     u = uM, copula = mC),
        398.66745, tol = 7e-7)
)

## "free" weights --- estimation == FIXME: will change after re-parametrization
optCtrl <- list(maxit = 1000, trace = TRUE)
if (doExtras) { # slowish
 try({
    st0 <- system.time(
        f. <- fitCopula(mC, uM, optim.method = "BFGS", optim.control=optCtrl, traceOpt=TRUE))
    ## fails: non-finite finite-difference value [7]
    print(st0) #
    print(lf. <- logLik(f.))
    print(summary(f.))
 })

    st1 <- system.time(
        ff <- fitCopula(mC, uM, optim.method = "Nelder-Mead", optim.control=optCtrl))
    ## warning: ... possible convergence problem: optim() gave code=1
    print(st1) # 11 sec
    print(lff <- logLik(ff))
    print(summary(ff))

    st2 <- system.time(
        f2 <- fitCopula(mC, uM, optim.method = "L-BFGS-B",    optim.control=optCtrl))
    print(st2) # 28 sec
    print(lf2 <- logLik(f2))
    print(summary(f2))
}

