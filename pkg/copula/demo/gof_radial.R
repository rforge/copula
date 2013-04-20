## Copyright (C) 2013 Marius Hofert
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


## Goal: Graphical goodness-of-fit test(s) for meta-elliptical and
##       meta-Archimedean models (or meta-'radial' models)


### setup ######################################################################

## load packages
require(Matrix)
require(mvtnorm)
require(copula)


### 1) Functions ###############################################################

##' @title Q-Q plots of angular distributions against Beta distributions
##' @param k k for which B_k is computed
##' @param Bmat matrix as returned by gofBTstat()
##' @return invisible()
##' @author Marius Hofert
qqp <- function(k, Bmat)
    qqplot2(Bmat[,k], qF=function(p) qbeta(p, shape1=k/2, shape2=(ncol(Bmat)+1-k)/2),
            main.args=list(text=as.expression(substitute(plain("Beta")(s1,s2)~~
                bold("Q-Q Plot"), list(s1=k/2, s2=(ncol(Bmat)+1-k)/2)))))

##' @title -log-likelihood for t copulas
##' @param nu d.o.f. parameter
##' @param P correlation matrix
##' @param u data matrix (in [0,1]^d)
##' @return -log-likelihood for a t copula
##' @author Marius Hofert
##' Note: requires 'mvtnorm'
nLLt <- function(nu, P, u){
    stopifnot((d <- ncol(u))==ncol(P), ncol(P)==nrow(P))
    qtu <- qt(u, df=nu)
    ldtnu <- function(u, P, nu) dmvt(qtu, sigma=P, df=nu, log=TRUE) -
        rowSums(dt(qtu, df=nu, log=TRUE)) # t copula log-density
    -sum(ldtnu(u, P=P, nu=nu))
}


### 2) Meta-elliptical models ##################################################

### 2.1) Generate data from a multivariate normal distribution #################

set.seed(100)
n <- 250
d <- 5

## build a mean vector and a covariance matrix, and generate multivariate normal data
mu <- d:1
L <- diag(d) # identity in dim d
L[lower.tri(L)] <- 1:(d*(d-1)/2)/d # Cholesky factor (diagonal entries have to be > 0)
Sigma <- L %*% t(L) # positive definite
X <- rmvnorm(n, mean=mu, sigma=Sigma) # multivariate normal data

## compute pseudo-observations of the radial part and uniform distribution
RUnorm <- RSpobs(X, method="ellip")
Rnorm <- RUnorm$R
Unorm <- RUnorm$U

## 2.1.1) Q-Q plot of R against the correct quantiles
qqplot2(Rnorm, qF=function(p) sqrt(qchisq(p, df=d)),
        main.args=list(text=as.expression(substitute(bold(italic(chi[d.])~~"Q-Q Plot"), list(d.=d))),
        side=3, cex=1.3, line=1.1, xpd=NA))

## 2.1.2) Q-Q plot of R against the quantiles of F_R for a multivariate t_4 distribution
nu <- 4
qqplot2(Rnorm, qF=function(p) sqrt(d*qf(p, df1=d, df2=nu)),
        main.args=list(text=as.expression(substitute(bold(italic(F[list(d.,nu.)](r^2/d.))~~"Q-Q Plot"),
            list(d.=d, nu.=nu))), side=3, cex=1.3, line=1.1, xpd=NA))


### 2.2) Generate data from a multivariate t4 distribution #####################

X <- rmvt(n, mean=mu, sigma=Sigma, df=nu) # multivariate t4 data
## note: cor(X) ~= cov2cor(Sigma), so Sigma is the *covariance* matrix of the t distribution

## compute pseudo-observations of the radial part and uniform distribution
RUe <- RSpobs(X, method="ellip")
Re <- RUe$R
Ue <- RUe$U

## 2.2.1) Q-Q plot of R against the correct quantiles
qqplot2(Re, qF=function(p) sqrt(d*qf(p, df1=d, df2=nu)),
        main.args=list(text=as.expression(substitute(bold(italic(F[list(d.,nu.)](r^2/d.))~~"Q-Q Plot"),
            list(d.=d, nu.=nu))), side=3, cex=1.3, line=1.1, xpd=NA))

## 2.2.2) Q-Q plot of R against the quantiles of F_R for a multivariate normal distribution
qqplot2(Re, qF=function(p) sqrt(qchisq(p, df=d)),
        main.args=list(text=as.expression(substitute(bold(italic(chi[d.])~~"Q-Q Plot"), list(d.=d))),
        side=3, cex=1.3, line=1.1, xpd=NA))

## 2.2.3) Q-Q plot of the angular distribution (Bmat[,k] should follow a Beta(k/2, (d-k)/2) distribution)
Bmat <- gofBTstat(Ue)
qqp(1, Bmat=Bmat) # k=1
qqp(3, Bmat=Bmat) # k=3

## 2.2.4) Check independence between radial part and B_1 and B_3
plot(pobs(cbind(Re, Bmat[,1])), xlab=expression(italic(R)), ylab=expression(italic(B)[1]),
     main=expression(bold("Rank plot between"~~italic(R)~~"and"~~italic(B)[1])))
plot(pobs(cbind(Re, Bmat[,3])), xlab=expression(italic(R)), ylab=expression(italic(B)[3]),
     main=expression(bold("Rank plot between"~~italic(R)~~"and"~~italic(B)[3])))


### 2.3) Generate data from a Clayton and a Gumbel copula ######################

tau <- 0.5
U.C <- rCopula(n, archmCopula("Clayton", param=getAcop("Clayton")@iTau(tau), dim=d))
U.G <- rCopula(n, archmCopula("Gumbel", param=getAcop("Gumbel")@iTau(tau), dim=d))

## compute pseudo-observations of the radial part
RUC <- RSpobs(U.C, method="ellip")
RC <- RUC$R
UC <- RUC$U
RUG <- RSpobs(U.G, method="ellip")
RG <- RUG$R
UG <- RUG$U

## 2.3.1) Q-Q plot of R against the quantiles of F_R for a multivariate t_4 distribution
qqplot2(RC, qF=function(p) sqrt(d*qf(p, df1=d, df2=nu)),
        main.args=list(text=as.expression(substitute(bold(italic(F[list(d.,nu.)](r^2/d.))~~"Q-Q Plot"),
            list(d.=d, nu.=nu))), side=3, cex=1.3, line=1.1, xpd=NA))

## 2.3.2) Q-Q plot of R against the quantiles of F_R for a multivariate t_4 distribution
qqplot2(RG, qF=function(p) sqrt(d*qf(p, df1=d, df2=nu)),
        main.args=list(text=as.expression(substitute(bold(italic(F[list(d.,nu.)](r^2/d.))~~"Q-Q Plot"),
            list(d.=d, nu.=nu))), side=3, cex=1.3, line=1.1, xpd=NA))

## 2.3.3) Q-Q plot of the angular distribution (Bmat[,k] should *not* follow a Beta(k/2, (d-k)/2) distribution anymore)
BmatC <- gofBTstat(UC)
BmatG <- gofBTstat(UG)
qqp(1, Bmat=BmatC) # k=1
qqp(3, Bmat=BmatC) # k=3
qqp(1, Bmat=BmatG) # k=1
qqp(3, Bmat=BmatG) # k=3

## 2.3.4) Check independence between radial part and B_1 and B_3
plot(pobs(cbind(RC, BmatC[,1])), xlab=expression(italic(R)), ylab=expression(italic(B)[1]),
     main=expression(bold("Rank plot between"~~italic(R)~~"and"~~italic(B)[1])))
plot(pobs(cbind(RG, BmatG[,1])), xlab=expression(italic(R)), ylab=expression(italic(B)[1]),
     main=expression(bold("Rank plot between"~~italic(R)~~"and"~~italic(B)[1])))


### 3) Meta-Archimedean ########################################################

set.seed(100)
n <- 250
d <- 5

## generate data from a meta-Gumbel model with N(0,1) margins
family <- "Gumbel"
tau <- 0.5
th <- iTau(archmCopula(family), tau)
copG <- archmCopula(family, param=th, dim=d)
U.G <- rCopula(n, copula=copG)
X <- qnorm(U.G)


RS <- RSpobs(X, method="archm")
R <- RS$R
S <- RS$S


### 4) Application #############################################################

data(SMI.12)
d <- ncol(x <- diff(log(SMI.12))) # log-returns
u <- pobs(x) # pseudo-observations
tau <- cor(u, method="kendall") # = cor(x, method="kendall")

plot TODO: visualize matrix of pairwise tau! => not Archimedean

## 4.1) estimate a multivariate t copula (with the approach of
##      Demarta, McNeil (2005)) assuming [and later checking] if
##      the data comes from a meta-t-model

## estimate P
P <- as.matrix(nearPD(sin(cor(x, method="kendall")*pi/2), corr=TRUE)$mat)

## estimate nu via MLE for given P
nu <- seq(.5, 64, by=.25)
nLL <- sapply(nu, nLLt, P=P, u=u)
plot(nu, nLL+1200, type="l", log="xy",
     main=expression(bold("Negative log-likelihood as a function in"~~nu)),
     xlab=bquote(nu), ylab=expression(1200-logL(nu)))
## now we got the picture, find the minimum:
nu. <- optimize(nLLt, interval=c(.5, 64), P=P, u=u, tol=1e-7)$minimum # 11.96

## 4.2) checking if the data indeed comes from a meta-t_{nu.}-model
RUe <- RSpobs(x, method="ellip")
Re <- RUe$R
Ue <- RUe$U

## 4.2.1) Q-Q plot of R against the quantiles of F_R for the estimated t distribution
qqplot2(Re, qF=function(p) sqrt(d*qf(p, df1=d, df2=nu.)),
        main.args=list(text=as.expression(substitute(bold(italic(F[list(d.,nu.)](r^2/d.))~~"Q-Q Plot"),
            list(d.=d, nu.=round(nu.,2)))), side=3, cex=1.3, line=1.1, xpd=NA))

## 4.2.2) Q-Q plot of the angular distribution (Bmat[,k] should follow a Beta(k/2, (d-k)/2) distribution)
Bmat <- gofBTstat(Ue)
qqp(1, Bmat=Bmat) # k=1
qqp(10, Bmat=Bmat) # k=10

## 4.2.3) Check independence between radial part and B_1 and B_3
plot(pobs(cbind(Re, Bmat[,1])), xlab=expression(italic(R)), ylab=expression(italic(B)[1]),
     main=expression(bold("Rank plot between"~~italic(R)~~"and"~~italic(B)[1])))
plot(pobs(cbind(Re, Bmat[,10])), xlab=expression(italic(R)), ylab=expression(italic(B)[10]),
     main=expression(bold("Rank plot between"~~italic(R)~~"and"~~italic(B)[10])))
