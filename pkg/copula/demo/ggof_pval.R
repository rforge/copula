## Copyright (C) 2012 Marius Hofert and Martin Maechler
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


### Pairwise bootstrapped  P-values  [expensively]:
### =================================================

require(copula)

## == Example 3: 5d t_4 copula  of demo(gof_graph) # ./gof_graph.R

if(!(exists("doX") && is.logical(as.logical(doX))))
    print(doX <- copula:::doExtras())
(N <- if(doX) 256 else 32)# be fast when run as "check"

## setup
n <- 1000 # sample size
d <- 5 # dimension
family <- "t" # copula family
df <- 4 # degrees of freedom
if(setSeeds) set.seed(4)

## define and sample the copula, build pseudo-observations
tau <- c(0.2, 0.4, 0.6)
r <- iTau(tCopula(), tau)
P <- c(r[2], r[1], r[1], r[1], # upper triangle (without diagonal) of correlation "matrix"
             r[1], r[1], r[1],
                   r[3], r[3],
                         r[3])
copt4 <- ellipCopula(family, param=P, dim=d, dispstr="un", df=df, df.fixed=TRUE)
U <- rCopula(n, cop=copt4)
U. <- pobs(U)

## define the H0 copula
## Note: that's the same result when using pseudo-observations since estimation via
##       tau is invariant strictly increasing transformations
stopifnot(require(Matrix))
P. <- nearPD(iTau(tCopula(), cor(U., method="kendall")))$mat # estimate P
P.. <- P2p(P.)
plot(P, P.., asp=1); abline(0,1, col=adjustcolor("gray", 0.9)) # P. should be close to P
copH0 <- ellipCopula(family, param=P.., dim=d, dispstr="un", df=df, df.fixed=TRUE)

## create array of pairwise copH0-transformed data columns
cu.u <- pairwiseCcop(U., copH0, df=df)

## compute pairwise matrix of p-values and corresponding colors
pwIT <- pairwiseIndepTest(cu.u, N=N, verbose=interactive()) # (d,d)-matrix of test results
round(pmat <- pviTest(pwIT), 3) # pick out p-values
cc <- pairsColList(pmat) # compute corresponding colors

## which pairs violate H0?
which(pmat < 0.05, arr.ind=TRUE) # [none]

## pairwise Rosenblatt plot
title <- list("Pairwise Rosenblatt transformed pseudo-observations",
              expression(bold("to test")~~italic(H[0]:C~~bold("is t")[4])))
main.line <- c(4, 1.4)
## --- JCGS, Fig.4(left) ---
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", main=title, main.line=main.line,
                main.centered=TRUE)
## --- JCGS, Fig.4(right) ---
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", main=title, main.line=main.line,
                method = "QQchisq", main.centered=TRUE)


##- Now compute the  P-values  very expensively, but
## "correctly" (well, if bootstrap is valid) :
