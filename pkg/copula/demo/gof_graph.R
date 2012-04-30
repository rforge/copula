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


### setup ######################################################################

require(copula)

source(system.file("Rsource", "wrapper.R", package="copula"))
source(system.file("Rsource", "gof-graph.R", package="copula"))
source(system.file("Rsource", "graphics.R", package="copula"))


### Example 1: 5d Gumbel copula ################################################

## setup
n <- 1000 # sample size
d <- 5 # dimension
family <- "Gumbel" # copula family

## define and sample the copula (= H0 copula), build pseudo-observations
cop <- getAcop(family)
th <- cop@tauInv(tau <- 0.5) # correct parameter value
copH0 <- onacopulaL(family, list(th, 1:d)) # define H0 copula
U <- rcop(n, cop=copH0)
U. <- pobs(U)

## create array of pairwise copH0-transformed data columns
cu.u <- pairwiseCcop(U., copH0)
stopifnot(is.array(cu.u), dim(cu.u) == c(n,d,d)) # check

## compute pairwise matrix of p-values and corresponding colors
pwIT <- pairwiseIndepTest(cu.u, verbose=TRUE) # (d,d)-matrix of test results
round(pmat <- pviTest(pwIT), 3) # pick out p-values
cc <- pairsColList(pmat) # compute corresponding colors

## which pairs violate H0?
which(pmat < 0.05, arr.ind=TRUE)

## check whether p-values are uniform
if(d > 10){
    n. <- d*(d-1)
    qqplot(qunif(ppoints(n.)), sort(pmat), main = paste("n = ",n.))
    abline(0,1)
}


### Example 1: Plots ###########################################################

## 1) plain
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".")

## 2) with title
pwRoto <- "Pairwise Rosenblatt transformed observations"
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", main=pwRoto)

## 3) with title and subtitle
gp <- format(gpviTest(pmat), digits=1, nsmall=1)
sub <- paste(names(gp), gp, sep=": ")
sub. <- paste(paste(sub[1:3], collapse=", "), "\n",
              paste(sub[4:6], collapse=", "), "\n",
              paste(sub[7:8], collapse=", "), sep="")
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", main=pwRoto, sub=sub., sub.line=6.8)

## 4) two-line title including expressions, and centered
title <- list(paste(pwRoto, "to test"),
              substitute(italic(H[0]:C~~bold("is Gumbel with"~~tau==tau.)),
                         list(tau.=tau)))
main.line <- c(4, 1.4)
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", main=title, main.line=main.line,
                main.centered=TRUE, sub=sub., sub.line=6.8)

## 5) omit panel borders
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", panel.border=FALSE)

## 6) without axes
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", axes=FALSE)

## 7) without axes and borders
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", axes=FALSE, panel.border=FALSE)

## 8) make black colors of the dots less dominant
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", col=adjustcolor("black", 0.5))

## 9) plot just colors (axis labels are automagically removed)
pairsRosenblatt(cu.u, pvalueMat=pmat, method="none")

## 10) also remove labels on the diagonal
pairsRosenblatt(cu.u, pvalueMat=pmat, method="none", labels="n")

## 11) Q-Q plots -- can, in general, better detect outliers
pairsRosenblatt(cu.u, pvalueMat=pmat, method="QQchisq", cex=0.2)
## pairsRosenblatt(cu.u, pvalueMat=pmat, method="QQchisq", cex=0.2,
##                 panel=function(x, y, ...){
##                     points(x, y, ...)
##                     qqline(y) ## only makes sense for the normal distribution
##                 })
## => TODO: maybe in the future with a more general function for qqline (supporting
##          other distributions)

## 12) P-P plots -- actually, MM sees *more* (though outliers are invisible)
pairsRosenblatt(cu.u, pvalueMat=pmat, method="PPchisq")
pairsRosenblatt(cu.u, pvalueMat=pmat, method="PPchisq",
                panel=function(x, y, ...){
                    points(x, y, ...)
                    abline(0, 1, col="green") # add straight line
                })


### Example 1: Boundary cases ##################################################

## Note: this is only for checking "boundary input cases", it does not make
##       sense given the data.

## 1) one pdiv, within range(pmat)
pmat <- matrix(c(
               NA, 0.1, 0.2, 0.3, 0.4,
               0.2, NA, 0.3, 0.4, 0.5,
               0.3, 0.4, NA, 0.5, 0.6,
               0.4, 0.5, 0.6, NA, 0.7,
               0.5, 0.6, 0.7, 0.8, NA
               ), nrow=5, ncol=5)
bgColList <- pairsColList(pmat, pdiv=0.2, signif.P=0.2,
                          colsBelow="green", colsAbove="orange")
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", bgColList=bgColList)

## 2) one pdiv, pdiv < min(pmat)
bgColList <- pairsColList(pmat, pdiv=0.05, signif.P=0.05,
                          colsBelow="green", colsAbove="orange") # note: we don't even need colsBelow here
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", bgColList=bgColList)

## 3) one pdiv, pdiv > max(pmat)
bgColList <- pairsColList(pmat, pdiv=0.9, signif.P=0.9, colsBelow="green", colsAbove="orange") # note: we don't even need colsAbove here
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", bgColList=bgColList)

## 4) one pdiv, equal to all values of pmat
pmat <- matrix(0.05, nrow=5, ncol=5)
diag(pmat) <- NA
bgColList <- pairsColList(pmat, pdiv=0.05, signif.P=0.05, colsBelow="green", colsAbove="orange") # note: we don't even need colsBelow here
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", bgColList=bgColList)
## => plot color key up to 1 (see pairsColList())

## 5) one p-value equal to 0; in this case we need pmin0 > 0
pmat[1,2] <- 0
bgColList <- pairsColList(pmat, signif.P=0.05, pmin0=1e-5)
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", bgColList=bgColList)

## 6) specifically call pairsColList but forget to set pmin0 to check the error message
if(FALSE){
    bgColList <- pairsColList(pmat, signif.P=0.05)
    pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", bgColList=bgColList)
}


### Example 2: (2,3)-nested Gumbel copula ######################################

## setup
n <- 1000 # sample size
d <- 5 # dimension
family <- "Gumbel" # copula family

## define and sample the copula, build pseudo-observations
cop <- getAcop(family)
th <- cop@tauInv(tau <- c(0.2, 0.4, 0.6))
nacList <- list(th[1], NULL, list(list(th[2], 1:2), list(th[3], 3:d)))
copG <- copCreate(family, nacList=nacList)
U <- rcop(n, cop=copG)
U. <- pobs(U)

## define the H0 copula
th0 <- cop@tauInv(tau0 <- c(0.2, 0.4, 0.4)) # wrong 2nd-sector-parameter
nacList <- list(th0[1], NULL, list(list(th0[2], 1:2), list(th0[3], 3:d)))
copH0 <- onacopulaL(family, nacList)

## create array of pairwise copH0-transformed data columns
cu.u <- pairwiseCcop(U., copH0)

## compute pairwise matrix of p-values and corresponding colors
pwIT <- pairwiseIndepTest(cu.u, verbose=TRUE) # (d,d)-matrix of test results
round(pmat <- pviTest(pwIT), 3) # pick out p-values
cc <- pairsColList(pmat) # compute corresponding colors

## which pairs violate H0?
which(pmat < 0.05, arr.ind=TRUE)

## pairwise Rosenblatt plot
title <- list(paste(pwRoto, "to test"),
              substitute(italic(H[0]:C~~bold("is nested Gumbel with"~~
                                                      tau[0]==tau0*","
                                                      ~~tau[1]==tau1*","
                                                      ~~tau[2]==tau2)),
                         list(tau0=tau0[1], tau1=tau0[2], tau2=tau0[3])))
main.line <- c(4, 1.4)
gp <- format(gpviTest(pmat), digits=1, nsmall=1)
sub <- paste(names(gp), gp, sep=": ")
sub. <- paste(paste(sub[1:3], collapse=", "), "\n",
              paste(sub[4:6], collapse=", "), "\n",
              paste(sub[7:8], collapse=", "), sep="")
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", main=title, main.line=main.line,
                sub=sub., sub.line=6.8)


### Example 3: 5d t_4 copula (fixed/known d.o.f., estimated P) #################

## setup
n <- 1000 # sample size
d <- 5 # dimension
family <- "t" # copula family
df <- 4 # degrees of freedom

## define and sample the copula, build pseudo-observations
tau <- c(0.2, 0.4, 0.6)
r <- tauInv(tau, family=family)
P <- c(r[2], r[1], r[1], r[1], # upper triangle (without diagonal) of correlation "matrix"
             r[1], r[1], r[1],
                   r[3], r[3],
                         r[3])
copt4 <- copCreate(family, theta=P, d=d, dispstr="un", df=df, df.fixed=TRUE)
U <- rcop(n, cop=copt4)
U. <- pobs(U)

## define the H0 copula
## Note: that's the same result when using pseudo-observations since estimation via
##       tau is invariant strictly increasing transformations
P. <- nearPD(tauInv(cor(U., method="kendall"), family=family))$mat # estimate P
P. <- P.[lower.tri(P.)] # note: upper.tri() would lead to the wrong result due to the required ordering
plot(P, P., asp=1); abline(0,1, col=adjustcolor("gray",0.5)) # P. should be close to P
copH0 <- copCreate(family, theta=P., d=d, dispstr="un", df=df, df.fixed=TRUE)

## create array of pairwise copH0-transformed data columns
cu.u <- pairwiseCcop(U., copH0, df=df)

## compute pairwise matrix of p-values and corresponding colors
pwIT <- pairwiseIndepTest(cu.u, verbose=TRUE) # (d,d)-matrix of test results
round(pmat <- pviTest(pwIT), 3) # pick out p-values
cc <- pairsColList(pmat) # compute corresponding colors

## which pairs violate H0?
which(pmat < 0.05, arr.ind=TRUE)

## pairwise Rosenblatt plot
title <- list("Pairwise Rosenblatt transformed pseudo-observations",
              expression(bold("to test")~~italic(H[0]:C~~bold("is t")[4])))
main.line <- c(4, 1.4)
gp <- format(gpviTest(pmat), digits=1, nsmall=1)
sub <- paste(names(gp), gp, sep=": ")
sub. <- paste(paste(sub[1:3], collapse=", "), "\n",
              paste(sub[4:6], collapse=", "), "\n",
              paste(sub[7:8], collapse=", "), sep="")
pairsRosenblatt(cu.u, pvalueMat=pmat, pch=".", main=title, main.line=main.line,
                sub=sub., sub.line=6.8)
