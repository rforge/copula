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

## basic settings
doPDF <- FALSE
doCrop <- TRUE


### 1) Functions ###############################################################

## start pdf() with nice defaults
start.pdf <- function(file="Rplots.pdf", doPDF=TRUE, width=6, height=6, ...)
{
    if(doPDF) pdf(file=file, width=width, height=height, ...)
    invisible()
}

## see the R package 'simsalapar' for a (much) more sophisticated version
dev.off.pdf <- function(file="Rplots.pdf", doPDF=TRUE, doCrop=TRUE, ...)
{
    dev.off(...)
    if(doCrop) {
        f <- file.path(getwd(), file)
        system(paste("pdfcrop --pdftexcmd pdftex", f, f, "1>/dev/null 2>&1"))
    }
    invisible()
}

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

### 2.1) Checking R ############################################################

### 2.1.1) Generate data from a multivariate normal / t distribution ###########

## setup
n <- 250 # sample size
d <- c(2, 10) # c(2, 10, 50, 100) # dimensions; numerically critical (due to estimator of P) for d > n/2
nu <- 1 # degrees of freedom
rho <- 0.5 # equi-correlation parameter
SigmaType <- c("EC", "AR") # equi-correlation, autoregressive

## go through the Sigma types and dimensions
for(st in SigmaType) {
    for(d. in d) {

        ## generate data
        set.seed(271) # set seed
        mu <- rep(0, d.) # mean (does not make a difference in our setup)
        Sigma <- switch(st,
                        "EC" = {
                            Sigma <- matrix(rep(rho, d.*d.), nrow=d., ncol=d.)
                            diag(Sigma) <- rep(1, d.)
                            Sigma
                        },
                        "AR" = {
                            outer(1:d., 1:d., FUN=function(i,j) rho^abs(i-j))
                        },
                        stop("wrong method"))
        X.norm <- rmvnorm(n, mean=mu, sigma=Sigma) # multivariate normal data
        X.t <- rep(mu, each=n) + rmvt(n, sigma=Sigma, df=nu) # multivariate t_nu data (Sigma = *dispersion* matrix)

        ## compute the (R, S) decomposition
        RS.norm.norm <- RSpobs(X.norm, method="ellip", qQg=qnorm)
        RS.norm.t <- RSpobs(X.norm, method="ellip", qQg=function(p) qt(p, df=nu))
        RS.t.norm <- RSpobs(X.t, method="ellip", qQg=qnorm)
        RS.t.t <- RSpobs(X.t, method="ellip", qQg=function(p) qt(p, df=nu))

        ## Q-Q plot: R normal against the correct quantiles
        file <- paste0("ggof_radial_true=normal_H0=normal_d=", d.,"_Sigma=", st, ".pdf")
        start.pdf(file=file, doPDF=doPDF)
        par(pty="s") # use a square plotting region
        qqplot2(RS.norm.norm$R, qF=function(p) sqrt(qchisq(p, df=d.)),
                main.args=list(text=as.expression(substitute(bold(italic(chi[d..])~~"Q-Q Plot"), list(d..=d.))),
                side=3, cex=1.3, line=1.1, xpd=NA))
        dev.off.pdf(file=file, doPDF=doPDF, doCrop=doCrop)

        ## Q-Q plot: R normal against the quantiles of F_g for a multivariate t_nu distribution
        file <- paste0("ggof_radial_true=normal_H0=t4_d=", d.,"_Sigma=", st, ".pdf")
        start.pdf(file=file, doPDF=doPDF)
        par(pty="s") # use a square plotting region
        qqplot2(RS.norm.t$R, qF=function(p) sqrt(d.*qf(p, df1=d., df2=nu)),
                main.args=list(text=as.expression(substitute(bold(italic(F[list(d..,nu.)](r^2/d..))~~"Q-Q Plot"),
                    list(d..=d., nu.=nu))), side=3, cex=1.3, line=1.1, xpd=NA))
        dev.off.pdf(file=file, doPDF=doPDF, doCrop=doCrop)

        ## Q-Q plot: R t_nu against the quantiles of F_g for a multivariate normal distribution
        file <- paste0("ggof_radial_true=t4_H0=normal_d=", d.,"_Sigma=", st, ".pdf")
        start.pdf(file=file, doPDF=doPDF)
        par(pty="s") # use a square plotting region
        qqplot2(RS.t.norm$R, qF=function(p) sqrt(qchisq(p, df=d.)),
                main.args=list(text=as.expression(substitute(bold(italic(chi[d..])~~"Q-Q Plot"), list(d..=d.))),
                side=3, cex=1.3, line=1.1, xpd=NA))
        dev.off.pdf(file=file, doPDF=doPDF, doCrop=doCrop)

        ## Q-Q plot: R t_nu against the correct quantiles
        file <- paste0("ggof_radial_true=t4_H0=t4_d=", d., "_Sigma=", st, ".pdf")
        start.pdf(file=file, doPDF=doPDF)
        par(pty="s") # use a square plotting region
        qqplot2(RS.t.t$R, qF=function(p) sqrt(d.*qf(p, df1=d., df2=nu)),
                main.args=list(text=as.expression(substitute(bold(italic(F[list(d..,nu.)](r^2/d..))~~"Q-Q Plot"),
                    list(d..=d., nu.=nu))), side=3, cex=1.3, line=1.1, xpd=NA))
        dev.off.pdf(file=file, doPDF=doPDF, doCrop=doCrop)

    }
}

## TODO: from here

require(mvtnorm)
require(copula)

d <- 10
rho <- 0.5
mu <- rep(0, d)
n <- 1000
Sigma <- outer(1:d, 1:d, FUN=function(i,j) rho^abs(i-j))
set.seed(271)
X.norm <- rmvnorm(n, mean=mu, sigma=Sigma) # multivariate normal data
X <- rep(mu, each=n) + rmvt(n, sigma=Sigma, df=10) # multivariate t data

RS.norm <- RSpobs(X, method="ellip", qQg=qnorm)
RS.t2  <- RSpobs(X, method="ellip", qQg=function(p) qt(p, df=2))
RS.t4  <- RSpobs(X, method="ellip", qQg=function(p) qt(p, df=4))
RS.t10 <- RSpobs(X, method="ellip", qQg=function(p) qt(p, df=10))

dens.norm <- density(RS.norm$R)
dens.t2 <- density(RS.t2$R)
dens.t4 <- density(RS.t4$R)
dens.t10 <- density(RS.t10$R)

xran <- range(dens.norm$x, dens.t2$x, dens.t4$x, dens.t10$x)
yran <- range(dens.norm$y, dens.t2$y, dens.t4$y, dens.t10$y)
## normal
plot(dens.norm, xlim=xran, ylim=yran, lty=2)
x <- dens.t2$x
lines(x, 2*x*dchisq(x^2, df=d))
## t2
lines(dens.t2, col="red", lty=2)
lines(x, 2*x*df(x^2/d, df1=d, df2=2)/d, col="red")
## t4
lines(dens.t4, col="darkgreen", lty=2) # estimate
lines(x, 2*x*df(x^2/d, df1=d, df2=4)/d, col="darkgreen") # true
## t10
lines(dens.t10, col="blue", lty=2) # estimate
lines(x, 2*x*df(x^2/d, df1=d, df2=10)/d, col="blue") # true



## ## 2.2.3) Q-Q plot of the angular distribution (Bmat[,k] should follow a Beta(k/2, (d-k)/2) distribution)
## Bmat <- gofBTstat(S.t)
## qqp(1, Bmat=Bmat) # k=1
## qqp(3, Bmat=Bmat) # k=3

## ## 2.2.4) Check independence between radial part and B_1 and B_3
## plot(pobs(cbind(R.t, Bmat[,1])), xlab=expression(italic(R)), ylab=expression(italic(B)[1]),
##      main=expression(bold("Rank plot between"~~italic(R)~~"and"~~italic(B)[1])))
## plot(pobs(cbind(R.t, Bmat[,3])), xlab=expression(italic(R)), ylab=expression(italic(B)[3]),
##      main=expression(bold("Rank plot between"~~italic(R)~~"and"~~italic(B)[3])))


### 2.1.2) Generate data from a Clayton, Frank, and Gumbel copula ##############

tau <- 0.5

## go through the Sigma types and dimensions
for(d. in d) {

    ## generate data
    set.seed(271) # set seed
    U.C <- rCopula(n, archmCopula("Clayton", param=getAcop("Clayton")@iTau(tau), dim=d.))
    U.F <- rCopula(n, archmCopula("Frank", param=getAcop("Frank")@iTau(tau), dim=d.))
    U.G <- rCopula(n, archmCopula("Gumbel", param=getAcop("Gumbel")@iTau(tau), dim=d.))

    ## compute the (R, S) decomposition
    RS.C.norm <- RSpobs(U.C, method="ellip", qQg=qnorm)
    RS.C.t <- RSpobs(U.C, method="ellip", qQg=function(p) qt(p, df=nu))
    RS.F.norm <- RSpobs(U.F, method="ellip", qQg=qnorm)
    RS.F.t <- RSpobs(U.F, method="ellip", qQg=function(p) qt(p, df=nu))
    RS.G.norm <- RSpobs(U.G, method="ellip", qQg=qnorm)
    RS.G.t <- RSpobs(U.G, method="ellip", qQg=function(p) qt(p, df=nu))

    ## Q-Q plot: R Clayton against the quantiles of F_g for a multivariate normal distribution
    file <- paste0("ggof_radial_true=C_H0=normal_d=", d.,".pdf")
    start.pdf(file=file, doPDF=doPDF)
    par(pty="s") # use a square plotting region
    qqplot2(RS.C.norm$R, qF=function(p) sqrt(qchisq(p, df=d.)),
            main.args=list(text=as.expression(substitute(bold(italic(chi[d..])~~"Q-Q Plot"), list(d..=d.))),
            side=3, cex=1.3, line=1.1, xpd=NA))
    dev.off.pdf(file=file, doPDF=doPDF, doCrop=doCrop)

    ## Q-Q plot: R Clayton against the quantiles of F_g for a multivariate t_nu distribution
    file <- paste0("ggof_radial_true=C_H0=t4_d=", d.,".pdf")
    start.pdf(file=file, doPDF=doPDF)
    par(pty="s") # use a square plotting region
    qqplot2(RS.C.t$R, qF=function(p) sqrt(d.*qf(p, df1=d., df2=nu)),
            main.args=list(text=as.expression(substitute(bold(italic(F[list(d..,nu.)](r^2/d..))~~"Q-Q Plot"),
                list(d..=d., nu.=nu))), side=3, cex=1.3, line=1.1, xpd=NA))
    dev.off.pdf(file=file, doPDF=doPDF, doCrop=doCrop)

    ## Q-Q plot: R Frank against the quantiles of F_g for a multivariate normal distribution
    file <- paste0("ggof_radial_true=F_H0=normal_d=", d.,".pdf")
    start.pdf(file=file, doPDF=doPDF)
    par(pty="s") # use a square plotting region
    qqplot2(RS.F.norm$R, qF=function(p) sqrt(qchisq(p, df=d.)),
            main.args=list(text=as.expression(substitute(bold(italic(chi[d..])~~"Q-Q Plot"), list(d..=d.))),
            side=3, cex=1.3, line=1.1, xpd=NA))
    dev.off.pdf(file=file, doPDF=doPDF, doCrop=doCrop)

    ## Q-Q plot: R Frank against the quantiles of F_g for a multivariate t_nu distribution
    file <- paste0("ggof_radial_true=F_H0=t4_d=", d.,".pdf")
    start.pdf(file=file, doPDF=doPDF)
    par(pty="s") # use a square plotting region
    qqplot2(RS.F.t$R, qF=function(p) sqrt(d.*qf(p, df1=d., df2=nu)),
            main.args=list(text=as.expression(substitute(bold(italic(F[list(d..,nu.)](r^2/d..))~~"Q-Q Plot"),
                list(d..=d., nu.=nu))), side=3, cex=1.3, line=1.1, xpd=NA))
    dev.off.pdf(file=file, doPDF=doPDF, doCrop=doCrop)

    ## Q-Q plot: R Gumbel against the quantiles of F_g for a multivariate normal distribution
    file <- paste0("ggof_radial_true=G_H0=normal_d=", d.,".pdf")
    start.pdf(file=file, doPDF=doPDF)
    par(pty="s") # use a square plotting region
    qqplot2(RS.G.norm$R, qF=function(p) sqrt(qchisq(p, df=d.)),
            main.args=list(text=as.expression(substitute(bold(italic(chi[d..])~~"Q-Q Plot"), list(d..=d.))),
            side=3, cex=1.3, line=1.1, xpd=NA))
    dev.off.pdf(file=file, doPDF=doPDF, doCrop=doCrop)

    ## Q-Q plot: R Gumbel against the quantiles of F_g for a multivariate t_nu distribution
    file <- paste0("ggof_radial_true=G_H0=t4_d=", d.,".pdf")
    start.pdf(file=file, doPDF=doPDF)
    par(pty="s") # use a square plotting region
    qqplot2(RS.G.t$R, qF=function(p) sqrt(d.*qf(p, df1=d., df2=nu)),
            main.args=list(text=as.expression(substitute(bold(italic(F[list(d..,nu.)](r^2/d..))~~"Q-Q Plot"),
                list(d..=d., nu.=nu))), side=3, cex=1.3, line=1.1, xpd=NA))
    dev.off.pdf(file=file, doPDF=doPDF, doCrop=doCrop)

}





## ## compute pseudo-observations of the radial part
## RS.C <- RSpobs(U.C, method="ellip")
## R.C <- RS.C$R
## S.C <- RS.C$S
## RS.G <- RSpobs(U.G, method="ellip")
## R.G <- RS.G$R
## S.G <- RS.G$S

## ## 2.3.1) Q-Q plot of R against the quantiles of F_R for a multivariate t_4 distribution
## qqplot2(R.C, qF=function(p) sqrt(d*qf(p, df1=d, df2=nu)),
##         main.args=list(text=as.expression(substitute(bold(italic(F[list(d.,nu.)](r^2/d.))~~"Q-Q Plot"),
##             list(d.=d, nu.=nu))), side=3, cex=1.3, line=1.1, xpd=NA))

## ## 2.3.2) Q-Q plot of R against the quantiles of F_R for a multivariate t_4 distribution
## qqplot2(R.G, qF=function(p) sqrt(d*qf(p, df1=d, df2=nu)),
##         main.args=list(text=as.expression(substitute(bold(italic(F[list(d.,nu.)](r^2/d.))~~"Q-Q Plot"),
##             list(d.=d, nu.=nu))), side=3, cex=1.3, line=1.1, xpd=NA))

## ## 2.3.3) Q-Q plot of the angular distribution (Bmat[,k] should *not* follow a Beta(k/2, (d-k)/2) distribution anymore)
## Bmat.C <- gofBTstat(S.C)
## Bmat.G <- gofBTstat(S.G)
## qqp(1, Bmat=Bmat.C) # k=1
## qqp(3, Bmat=Bmat.C) # k=3
## qqp(1, Bmat=Bmat.G) # k=1
## qqp(3, Bmat=Bmat.G) # k=3

## ## 2.3.4) Check independence between radial part and B_1 and B_3
## plot(pobs(cbind(R.C, Bmat.C[,1])), xlab=expression(italic(R)), ylab=expression(italic(B)[1]),
##      main=expression(bold("Rank plot between"~~italic(R)~~"and"~~italic(B)[1])))
## plot(pobs(cbind(R.G, Bmat.G[,1])), xlab=expression(italic(R)), ylab=expression(italic(B)[1]),
##      main=expression(bold("Rank plot between"~~italic(R)~~"and"~~italic(B)[1])))


### 3) Meta-Archimedean ########################################################

if(FALSE) {
    ## checking qacR()
    set.seed(271)
    n <- 250
    d <- 5
    th <- 2
    family <- "Gumbel"
    p <- ppoints(n)
    qR <- qacR(p, family=family, theta=th, d=d, interval=c(0, 200)) # ~ 20s
    p. <- pacR(qR, family=family, theta=th, d=d) # check
    summary(p-p.) # => fine
}

## generate data from a meta-Gumbel model with N(0,1) margins
family <- "Gumbel"
tau <- 0.5
th <- iTau(archmCopula(family), tau)
copG <- archmCopula(family, param=th, dim=d)
set.seed(271)
U.G <- rCopula(n, copula=copG)
X <- qnorm(U.G)

## (R, S) decomposition
RS.G <- RSpobs(X, method="archm")

## Q-Q plot: R Gumbel against the correct quantiles
file <- paste0("ggof_radial_true=G_2_H0=.pdf")
start.pdf(file=file, doPDF=doPDF)
par(pty="s") # use a square plotting region
qqplot2(RS.G$R, qF=function(p) qacR(p, family=family, theta=th, d=d, interval=c(0, 1e12)),
        main.args=list(text=as.expression(substitute(bold(italic(F[list(d..,nu.)](r^2/d..))~~"Q-Q Plot"),
            list(d..=d., nu.=nu))), side=3, cex=1.3, line=1.1, xpd=NA))
dev.off.pdf(file=file, doPDF=doPDF, doCrop=doCrop)


### 4) Application #############################################################

data(SMI.12)
d <- ncol(x <- diff(log(SMI.12))) # log-returns
u <- pobs(x) # pseudo-observations
tau <- cor(u, method="kendall") # = cor(x, method="kendall")

## plot TODO: visualize matrix of pairwise tau! => not Archimedean

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
RS.t <- RSpobs(x, method="ellip")
R.t <- RS.t$R
S.t <- RS.t$S

## 4.2.1) Q-Q plot of R against the quantiles of F_R for the estimated t distribution
qqplot2(R.t, qF=function(p) sqrt(d*qf(p, df1=d, df2=nu.)),
        main.args=list(text=as.expression(substitute(bold(italic(F[list(d.,nu.)](r^2/d.))~~"Q-Q Plot"),
            list(d.=d, nu.=round(nu.,2)))), side=3, cex=1.3, line=1.1, xpd=NA))

## 4.2.2) Q-Q plot of the angular distribution (Bmat[,k] should follow a Beta(k/2, (d-k)/2) distribution)
Bmat <- gofBTstat(S.t)
qqp(1, Bmat=Bmat) # k=1
qqp(10, Bmat=Bmat) # k=10

## 4.2.3) Check independence between radial part and B_1 and B_3
plot(pobs(cbind(R.t, Bmat[,1])), xlab=expression(italic(R)), ylab=expression(italic(B)[1]),
     main=expression(bold("Rank plot between"~~italic(R)~~"and"~~italic(B)[1])))
plot(pobs(cbind(R.t, Bmat[,10])), xlab=expression(italic(R)), ylab=expression(italic(B)[10]),
     main=expression(bold("Rank plot between"~~italic(R)~~"and"~~italic(B)[10])))
