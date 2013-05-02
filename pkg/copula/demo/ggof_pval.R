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

if(!(exists("setSeeds") && is.logical(as.logical(setSeeds))))
setSeeds <- TRUE # for reproducibility
##  maybe set to FALSE *before* running this demo

if(!(exists("doX") && is.logical(as.logical(doX))))
    print(doX <- copula:::doExtras())
(N <- if(doX) 256 else 32)# be fast when run as "check"

## setup
n <- 1000 # sample size
d <- 5 # dimension
family <- "t" # copula family
df <- 4 # degrees of freedom

## define and sample the copula, build pseudo-observations
tau <- c(0.2, 0.4, 0.6)
r <- iTau(tCopula(), tau)
P <- c(r[2], r[1], r[1], r[1], # upper triangle (without diagonal) of correlation "matrix"
             r[1], r[1], r[1],
                   r[3], r[3],
                         r[3])
copt4 <- ellipCopula(family, param=P, dim=d, dispstr="un", df=df, df.fixed=TRUE)
if(setSeeds) set.seed(4)
U <- rCopula(n, cop=copt4)
U. <- pobs(U)

stopifnot(require(Matrix))

## define the H0 copula
## Note: that's the same result when using pseudo-observations since estimation via
##       tau is invariant strictly increasing transformations
mkH0cop <- function(u, df) {
    stopifnot((d <- ncol(u)) >= 2)
    P. <- nearPD(iTau(tCopula(), cor(u, method="kendall")))$mat # estimate P
    ellipCopula("t", param= P2p(P.), dim = d, dispstr="un", df=df, df.fixed=TRUE)
}

copH0 <- mkH0cop(U., df=df)

iTest <- indepTestSim(n, p=2, m=2, N=N, verbose = TRUE)

##' Array of pairwise  copH0-transformed "data", and P-values
mkcP <- function(u, df, copH0 = mkH0cop(u, df=df), N, verbose=interactive(),
                 P.only=FALSE)
{
    ## create array of pairwise copH0-transformed data columns
    cu.u <- pairwiseCcop(u, copH0, df=df)
    ## compute pairwise matrix of p-values and corresponding colors
    pwIT <- pairwiseIndepTest(cu.u, N=N, verbose=verbose) # (d,d)-matrix of test results
    ## pick out p-values:
    pwIT <- pviTest(pwIT)
    if(P.only) pwIT else list(cu.u=cu.u, P = pwIT)
}

## compute pairwise matrix of p-values and corresponding colors
cP <- mkcP(U., df=df, copH0=copH0, N=N) # pick out p-values

round(pmat <- cP$P, 3)
cc <- pairsColList(pmat) # compute corresponding colors
cu.u <- cP[["cu.u"]]

## which pairs violate H0?
which(pmat < 0.05, arr.ind=TRUE) # [none]

## pairwise Rosenblatt plot
title <- list("Pairwise Rosenblatt transformed pseudo-observations",
              expression(bold("to test")~~italic(H[0]:C~~bold("is t")[4])))
line.main <- c(4, 1.4)
## --- JCGS, Fig.4(left) ---
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", main=title, line.main=line.main,
                main.centered=TRUE)
## --- JCGS, Fig.4(right) ---
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", main=title, line.main=line.main,
                method = "QQchisq", main.centered=TRUE)


##- Now compute the  P-values  very expensively, but
## "correctly" (well, if bootstrap is valid) :

(N <- if(doX) 128 else 16)# of "internal" bootstrap
B1 <- 500
B2 <- 256

## For testing, just during developement, take very small numbers:  [FIXME]
N <-  8; op <- options("copula:warn.idTS" = FALSE)
N <- 32; op <- options("copula:warn.idTS" = FALSE)
B1 <- 3
B2 <- 5

trace.lev <- 2

PP <- array(NA_real_, dim= c(d,d, B1, B2))
if(setSeeds)
set.seed(17)
if(trace.lev >= 1) cat(sprintf("OUTER bootstrap (b in 1:%d): \n", B1))
for(b in 1:B1) {
    ii <- sample.int(n, replace=TRUE)
    Us <- U.[ii, ] ## U*(b) {pseudo obs}
    ## estimate theta_n  and  H0_{s, theta_n}  (via iTau):
    copH0 <- mkH0cop(Us, df=df)
    if(trace.lev >= 1) cat(sprintf("b=%d, Det(\\hat{P})= %g\n", b,
       det(p2P(copH0@parameters, d = copH0@dimension))))
    ##
    ##
    ## Now for this  simple null hypothesis, need to boostrap *again*:
    if(trace.lev >= 2) cat(sprintf(" inner bootstrap (jj in 1:%d): ", B2))
    for(jj in 1:B2) {
        if(trace.lev >= 2) cat(if(jj %% 50 == 0) ".\n    " else
           if(jj %% 10 == 0) (jj / 10) %% 10 else ".")
        i2 <- sample.int(n, replace=TRUE)
        U2 <- Us[i2, ]
        ##
        ##
        ## ....
        ## ....
        ##
        ## and get P-values for it:
	PP[ , , b, jj] <- mkcP(U2, df=df, copH0=copH0, N=N, P.only=TRUE,
                               verbose = (trace.lev >= 3))
        ##
    }## for(jj in 1:B2)  inner bootstrap
    if(trace.lev >= 2) cat("\n")
    ##
    ## ...
    ## ...
    ##
}## for(b in 1:B1) --- *outer* bootstrap

str(PP)
## FIXME?  only store *means* of P values;
## FIXME?  or   store even the parameter vectors in addition?

options(op)# revert
