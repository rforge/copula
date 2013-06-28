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


## Goal: Graphical goodness-of-fit test(s) for meta-Archimedean and
##       meta-elliptical models (or meta-'radial' models)


## TODO:
## - shall we use pobs() everywhere? so far: known margins
## - 1.2), 1.3): map to chi^2 (1d setup)?


### Setup ######################################################################

## load packages and sources
require(Matrix)
require(mvtnorm)
require(copula)
source(system.file("Rsource", "AC-Liouville.R", package="copula"))

## basic settings
.seed <- 271 # for seeding
doPDF <- TRUE # whether plotting is done to pdf; TODO FALSE



### 0) Minor checks ############################################################

if(FALSE) {

    ## checking pacR() (F_R df) for Clayton
    family <- "Clayton"
    tau <- 0.5
    m <- 256
    dmax <- 20
    x <- seq(0, 20, length.out=m)
    y <- vapply(1:dmax, function(d)
                pacR(x, family=family, theta=iTau(archmCopula(family), tau), d=d),
                rep(NA_real_, m))
    plot(x, y[,1], type="l", ylim=c(0,1),
         xlab = expression(italic(x)),
         ylab = substitute(italic(F[R](x))~~"for d=1:"*dm, list(dm=dmax)))
    for(k in 2:dmax) lines(x, y[,k])

    ## checking qacR()
    set.seed(.seed)
    n <- 250
    d <- 5
    th <- 2
    family <- "Gumbel"
    p <- ppoints(n)
    qR <- qacR(p, family=family, theta=th, d=d, interval=c(0, 200)) # ~ 20s
    p. <- pacR(qR, family=family, theta=th, d=d) # check
    summary(p-p.) # => fine

}



### 1) Archimedean models ######################################################

## Data from a Clayton, Gumbel, and Gamma-Archimedean (R ~ Gamma) copula

## auxiliary functions

## Q-Q plot for R against the F_R quantiles for Clayton
qq.R.C <- function(x, d, tau, doPDF, file)
    qqplot2(x, qF = function(p) qacR(p, family="Clayton",
               theta = iTau(claytonCopula(), tau=tau),
               d=d, interval = c(1e-4, 1e2)), log="xy",
            main.args=list(text=expression(bold(italic(F[R]^(-1))~~"Q-Q Plot"~~
                "for"~~italic(F[R])~~"from a Clayton copula")),
            side=3, cex=1.3, line=1.1, xpd=NA), doPDF=doPDF, file=file)

## Q-Q plot for R against the F_R quantiles for Gumbel
qq.R.G <- function(x, d, tau, doPDF, file)
    qqplot2(x, qF = function(p) qacR(p, family="Gumbel",
               theta = iTau(gumbelCopula(), tau=tau),
               d=d, interval = c(1+1e-4, 1e2)), log="xy",
            main.args=list(text=expression(bold(italic(F[R]^(-1))~~"Q-Q Plot"~~
                "for"~~italic(F[R])~~"from a Gumbel copula")),
            side=3, cex=1.3, line=1.1, xpd=NA), doPDF=doPDF, file=file)

## Q-Q plot for R against Gamma quantiles
qq.R.g <- function(x, d, tau, doPDF, file)
{
    th <- iTauACsimplex(tau, d=d, Rdist="Gamma", interval=c(1e-2, 1e2))
    qqplot2(x, qF = function(p) qgamma(p, shape=th), log="xy",
            main.args=list(text=expression(bold(italic(F[R]^(-1))~~"Q-Q Plot"~~
                "for"~~italic(F[R])~~"being Gamma")),
            side=3, cex=1.3, line=1.1, xpd=NA), doPDF=doPDF, file=file)
}

## Q-Q plot for S
qq.angular <- function(x, shape1, shape2, doPDF, file)
    qqplot2(x, qF=function(p) qbeta(p, shape1=1, shape2=d),
            main.args=list(text=as.expression(substitute(plain("Beta")(s1,s2)~~
                bold("Q-Q Plot"), list(s1=1, s2=d)))),
            doPDF=doPDF, file=file)

## check independence of R and S (via B)
indepRB <- function(R, B, index, doPDF, file, width=6, height=6, crop=NULL, ...)
{
    if(doPDF) pdf(file=file, width=width, height=height)
    ## plot
    plot(pobs(cbind(R, B)),
         xlab=expression(italic(R)), ylab=substitute(italic(B)[i], list(i=index)),
         main=substitute(bold("Rank plot between"~~italic(R)~~"and"~~italic(B)[i]),
         list(i=index)), ...)
    ## pdf device
    if(doPDF) {
        r <- dev.off(...)
        iNcrop <- is.null(crop)
        if(.Platform$OS.type != "unix" && iNcrop) {
            warning("'crop = NULL' is only suitable for Unix") # => continue without cropping
        } else { # cropping
            f <- file.path(getwd(), file)
            if(iNcrop) { # crop with default command
                system(paste("pdfcrop --pdftexcmd pdftex", f, f, "1>/dev/null 2>&1"))
            } else if(nzchar(crop)) { # crop != "" crop with provided command
                system(crop)
            }
        }
        invisible(r)
    } else invisible() # return invisibly
}

## check independence of R and S (via mapping R to Gamma(d))
indepRS <- function(R, S, doPDF, file, width=6, height=6, crop=NULL, ...)
{
    d <- ncol(X)
    qR <- qgamma(rank(R)/(n+1), shape=d) # make it Gamma(d) (=> C=Pi)
    if(doPDF) pdf(file=file, width=width, height=height)
    ## plot
    pairs(pobs(qR*S), gap=0,
          labels=as.expression( sapply(seq_len(d), function(j) bquote(italic(tilde(R)S[.(j)]))) ),
          ...)
    ## pdf device
    if(doPDF) {
        r <- dev.off(...)
        iNcrop <- is.null(crop)
        if(.Platform$OS.type != "unix" && iNcrop) {
            warning("'crop = NULL' is only suitable for Unix") # => continue without cropping
        } else { # cropping
            f <- file.path(getwd(), file)
            if(iNcrop) { # crop with default command
                system(paste("pdfcrop --pdftexcmd pdftex", f, f, "1>/dev/null 2>&1"))
            } else if(nzchar(crop)) { # crop != "" crop with provided command
                system(crop)
            }
        }
        invisible(r)
    } else invisible() # return invisibly
}

## htrafo() for (R,S) data
htrafoRS <- function(R,S)
{
    ## essentially the same as in 'copula'
    lpsiI <- log(S)
    lcumsum <- matrix(unlist(lapply(1:d, function(j)
                                    copula:::lsum(t(lpsiI[,1:j, drop=FALSE])))),
                      ncol=d)
    U. <- matrix(unlist(lapply(1:(d-1),
                               function(k) exp(k*(lcumsum[,k]-
                                                  lcumsum[,k+1])) )),
                 ncol=d-1)
    cbind(U., rank(R)/(n+1))
}

## pairs plot for Hering--Hofert transformed data
pairs2 <- function(x, gap=0, doPDF, file, width=6, height=6, crop=NULL, ...)
{
    ## plot
    pairs(x, gap=0,
          labels=as.expression( sapply(seq_len(d), function(j) bquote(italic(U*"'"[.(j)]))) ),
          ...)
    ## pdf device
    if(doPDF) {
        r <- dev.off(...)
        iNcrop <- is.null(crop)
        if(.Platform$OS.type != "unix" && iNcrop) {
            warning("'crop = NULL' is only suitable for Unix") # => continue without cropping
        } else { # cropping
            f <- file.path(getwd(), file)
            if(iNcrop) { # crop with default command
                system(paste("pdfcrop --pdftexcmd pdftex", f, f, "1>/dev/null 2>&1"))
            } else if(nzchar(crop)) { # crop != "" crop with provided command
                system(crop)
            }
        }
        invisible(r)
    } else invisible() # return invisibly
}


## main
ggofArch <- function(n, d, tau, doPDF)
{
    ## generate data
    set.seed(.seed) # set seed
    U.C <- rCopula(n, archmCopula("Clayton", param = getAcop("Clayton")@iTau(tau),
                                  dim = d))
    U.G <- rCopula(n, archmCopula("Gumbel", param = getAcop("Gumbel")@iTau(tau),
                                  dim = d))
    U.g <- rACsimplex(n, d=d, theta = iTauACsimplex(tau, d=d, Rdist="Gamma",
                                                    interval=c(1e-2, 1e2))
                      Rdist="Gamma")

    ## compute the (R,S) decomposition for all data sets
    RS.C <- RSpobs(U.C, method="archm")
    RS.G <- RSpobs(U.G, method="archm")
    RS.g <- RSpobs(U.g, method="archm")

    ## plot options
    par(pty="s")


    ## 1.1) R-S decomposition ##################################################

    ## 1.1.1) Checking the radial part R #######################################

    ## data: R from Clayton...
    ## ... Q-Q plot against the F_R quantiles for Clayton (correct quantiles)
    file <- paste0("ggof_R_H0=C_true=C_d=", d, "_tau=", tau, ".pdf")
    qq.R.C(RS.C$R, d=d, tau=tau, doPDF=doPDF, file=file)
    ## ... Q-Q plot against the F_R quantiles for Gumbel
    file <- paste0("ggof_R_H0=G_true=C_d=", d, "_tau=", tau, ".pdf")
    qq.R.G(RS.C$R, d=d, tau=tau, doPDF=doPDF, file=file)
    ## ... Q-Q plot against Gamma quantiles
    file <- paste0("ggof_R_H0=Gamma_true=C_d=", d, "_tau=", tau, ".pdf")
    qq.R.g(RS.C$R, d=d, tau=tau, doPDF=doPDF, file=file)

    ## data: R from Gumbel...
    ## ... Q-Q plot against the F_R quantiles for Clayton
    file <- paste0("ggof_R_H0=C_true=G_d=", d, "_tau=", tau, ".pdf")
    qq.R.C(RS.G$R, d=d, tau=tau, doPDF=doPDF, file=file)
    ## ... Q-Q plot against the F_R quantiles for Gumbel (correct quantiles)
    file <- paste0("ggof_R_H0=G_true=G_d=", d, "_tau=", tau, ".pdf")
    qq.R.G(RS.G$R, d=d, tau=tau, doPDF=doPDF, file=file)
    ## ... Q-Q plot against Gamma quantiles
    file <- paste0("ggof_R_H0=Gamma_true=G_d=", d, "_tau=", tau, ".pdf")
    qq.R.g(RS.G$R, d=d, tau=tau, doPDF=doPDF, file=file)

    ## data: R ~ Gamma...
    ## ... Q-Q plot against the F_R quantiles for Clayton
    file <- paste0("ggof_R_H0=C_true=Gamma_d=", d, "_tau=", tau, ".pdf")
    qq.R.C(RS.g$R, d=d, tau=tau, doPDF=doPDF, file=file)
    ## ... Q-Q plot against the F_R quantiles for Gumbel
    file <- paste0("ggof_R_H0=G_true=Gamma_d=", d, "_tau=", tau, ".pdf")
    qq.R.G(RS.g$R, d=d, tau=tau, doPDF=doPDF, file=file)
    ## ... Q-Q plot against Gamma quantiles (correct quantiles)
    file <- paste0("ggof_R_H0=Gamma_true=Gamma_d=", d, "_tau=", tau, ".pdf")
    qq.R.g(RS.g$R, d=d, tau=tau, doPDF=doPDF, file=file)


    ## 1.1.2) Checking the angular part S ######################################

    ## Note (see Devroye (1986, pp. 207)):
    ##      B_j = S_1+..+S_j = U_{(j)} ~ Beta(j, d-j+1)
    ## We use j=1 and, if d>2, j=floor(d/2)
    d2 <- floor(d/2)

    ## data: S from Clayton
    ## B_1 (= S_1)
    file <- paste0("ggof_B1_true=C_d=", d, "_tau=", tau, ".pdf")
    qq.angular(RS.C$S[,1], shape1=1, shape2=d, doPDF=doPDF, file=file)
    ## B_{d/2}
    if(d > 2) {
        Sd2.C <- rowSums(RS.C$S[,seq_len(d2)])/rowSums(RS.C$S)
        file <- paste0("ggof_B", d2, "_true=C_d=", d, "_tau=", tau, ".pdf")
        qq.angular(Sd2.C, shape1=d2, shape2=d-d2+1, doPDF=doPDF, file=file)
    }

    ## data: S from Gumbel
    ## B_1 (= S_1)
    file <- paste0("ggof_B1_true=G_d=", d, "_tau=", tau, ".pdf")
    qq.angular(RS.G$S[,1], shape1=1, shape2=d, doPDF=doPDF, file=file)
    ## B_{d/2}
    if(d > 2) {
        Sd2.G <- rowSums(RS.G$S[,seq_len(d2)])/rowSums(RS.G$S)
        file <- paste0("ggof_B", d2, "_true=C_d=", d, "_tau=", tau, ".pdf")
        qq.angular(Sd2.G, shape1=d2, shape2=d-d2+1, doPDF=doPDF, file=file)
    }

    ## data: S from Gamma-R Archimedean copula
    ## B_1 (= S_1)
    file <- paste0("ggof_B1_true=Gamma_d=", d, "_tau=", tau, ".pdf")
    qq.angular(RS.g$S[,1], shape1=1, shape2=d, doPDF=doPDF, file=file)
    ## B_{d/2}
    if(d > 2) {
        Sd2.g <- rowSums(RS.g$S[,seq_len(d2)])/rowSums(RS.g$S)
        file <- paste0("ggof_B", d2, "_true=Gamma_d=", d, "_tau=", tau, ".pdf")
        qq.angular(Sd2.g, shape1=d2, shape2=d-d2+1, doPDF=doPDF, file=file)
    }


    ## 1.1.3) Checking independence of R and B_1, B_{d/2} ######################

    ## data: S from Clayton
    file <- paste0("ggof_indep_R_B1_true=C_d=", d, "_tau=", tau, ".pdf")
    indepRB(RS.C$R, B=RS.C$S[,1], index=1, doPDF=doPDF, file=file)
    file <- paste0("ggof_indep_R_B", d2, "_true=C_d=", d, "_tau=", tau, ".pdf")
    indepRB(RS.C$R, B=Sd2.C, index=d2, doPDF=doPDF, file=file)

    ## data: S from Gumbel
    file <- paste0("ggof_indep_R_B1_true=G_d=", d, "_tau=", tau, ".pdf")
    indepRB(RS.G$R, B=RS.G$S[,1], index=1, doPDF=doPDF, file=file)
    file <- paste0("ggof_indep_R_B", d2, "_true=C_d=", d, "_tau=", tau, ".pdf")
    indepRB(RS.G$R, B=Sd2.G, index=d2, doPDF=doPDF, file=file)

    ## data: S from Gamma-R Archimedean copula
    file <- paste0("ggof_indep_R_B1_true=Gamma_d=", d, "_tau=", tau, ".pdf")
    indepRB(RS.g$R, B=RS.G$S[,1], index=1, doPDF=doPDF, file=file)
    file <- paste0("ggof_indep_R_B", d2, "_true=C_d=", d, "_tau=", tau, ".pdf")
    indepRB(RS.g$R, B=Sd2.G, index=d2, doPDF=doPDF, file=file)


    ## 1.2) Other idea for checking independence between R and S ###############

    ## data: R from Clayton
    file <- paste0("ggof_indep_R_S_true=C_d=", d, "_tau=", tau, ".pdf")
    indepRS(RS.C$R, S=RS.C$S, doPDF=doPDF, file=file, pch=".") # should be independent

    ## data: R from Gumbel
    file <- paste0("ggof_indep_R_S_true=G_d=", d, "_tau=", tau, ".pdf")
    indepRS(RS.G$R, S=RS.G$S, doPDF=doPDF, file=file, pch=".") # should be independent

    ## data: R ~ Gamma
    file <- paste0("ggof_indep_R_S_true=Gamma_d=", d, "_tau=", tau, ".pdf")
    indepRS(RS.g$R, S=RS.g$S, doPDF=doPDF, file=file, pch=".") # should be independent


    ## 1.3) Hofert-Hering transform as GoF test ################################

    ## data: R from Clayton
    file <- paste0("ggof_htrafo_true=C_d=", d, "_tau=", tau, ".pdf")
    pairs2(htrafoRS(RS.C$R, S=RS.C$S), gap=0, doPDF=doPDF, file=file)

    ## data: R from Gumbel
    file <- paste0("ggof_htrafo_true=G_d=", d, "_tau=", tau, ".pdf")
    pairs2(htrafoRS(RS.G$R, S=RS.G$S), gap=0, doPDF=doPDF, file=file)

    ## data: R ~ Gamma
    file <- paste0("ggof_htrafo_true=Gamma_d=", d, "_tau=", tau, ".pdf")
    pairs2(htrafoRS(RS.g$R, S=RS.g$S), gap=0, doPDF=doPDF, file=file)
}

## setup
n <- 250 # sample size
d <- c(2, 10, 50) # dimensions
tau <- 0.5 # Kendall's tau

## call
ggofArch(n, d=d[1], tau=tau, doPDF=doPDF)
ggofArch(n, d=d[2], tau=tau, doPDF=doPDF)
ggofArch(n, d=d[3], tau=tau, doPDF=doPDF)



### 2) Exchangeable (but not Archimedean) models ###############################

## main
ggofExch <- function(n, d, tau, doPDF)
{
    ## generate data
    set.seed(.seed) # set seed
    rho <- iTau(normalCopula(), tau=tau)
    U.Ga <- rCopula(n, ellipCopula("normal", param = rho, dim = d))
    U.t4 <- rCopula(n, ellipCopula("t", param = rho, dim = d, df = 4))
    th.gL <- if(tau==0.5) 1.5 else # roughly 1.5 for tau==0.5 (we make it 'fast' here)
             iTauLiouville(tau, alpha=c(4,4), interval=c(1, 10), n.MC=1e5)
    U.gL <- rLiouville(n, alpha=rep(4, d), theta = th.gL, Rdist="Gamma")
    if(FALSE)
        pairs(U.gL, gap=0, pch=".") # check

    ## compute the (R,S) decomposition for all data sets
    RS.Ga <- RSpobs(U.Ga, method="archm")
    RS.t4 <- RSpobs(U.t4, method="archm")
    RS.gL <- RSpobs(U.gL, method="archm")

    ## plot options
    par(pty="s")

    ## TODO

}

## setup
n <- 250 # sample size
d <- c(2, 10, 50) # dimensions
tau <- 0.5 # Kendall's tau

## call
ggofExch(n, d=d[1], tau=tau, doPDF=doPDF)
ggofExch(n, d=d[2], tau=tau, doPDF=doPDF)
ggofExch(n, d=d[3], tau=tau, doPDF=doPDF)

2.2) check R with RS-decomposition
...
2.3) map hering hofert
Was meint johanna mit unterer 2. Tafel?



## setup
n <- 250 # sample size
d <- c(2, 10) # c(2, 10, 50, 100) # dimensions; numerically critical (due to estimator of P) for d > n/2
nu <- 1 # degrees of freedom
rho <- 0.5 # equi-correlation parameter
SigmaType <- c("EC", "AR") # equi-correlation, autoregressive

## go through the Sigma types and dimensions
for(st in SigmaType) {
    for(d. in d) {

        ## generate multivariate normal and t data
        set.seed(.seed) # set seed
        mu <- rep(0, d.) # mean (does not make a difference in our setup)
        Sigma <- switch(st, # standardized dispersion matrix
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
        X.t <- rep(mu, each=n) + rmvt(n, sigma=Sigma, df=nu) # multivariate t_nu data

        ## compute the (R,S) decomposition
        RS.norm.norm <- RSpobs(X.norm, method="ellip", qQg=qnorm)
        RS.norm.t <- RSpobs(X.norm, method="ellip", qQg=function(p) qt(p, df=nu))
        RS.t.norm <- RSpobs(X.t, method="ellip", qQg=qnorm)
        RS.t.t <- RSpobs(X.t, method="ellip", qQg=function(p) qt(p, df=nu))

        ## plot options
        par(pty="s")

        ## Q-Q plot: R normal against the correct quantiles
        file <- paste0("ggof_radial_true=normal_H0=normal_d=", d.,"_Sigma=", st, ".pdf")
        qqplot2(RS.norm.norm$R, qF=function(p) sqrt(qchisq(p, df=d.)),
                main.args=list(text=as.expression(substitute(bold(italic(chi[d..])~~"Q-Q Plot"), list(d..=d.))),
                side=3, cex=1.3, line=1.1, xpd=NA),
                doPDF=doPDF, file=file)

        ## Q-Q plot: R normal against the quantiles of F_g for a multivariate t_nu distribution
        file <- paste0("ggof_radial_true=normal_H0=t4_d=", d.,"_Sigma=", st, ".pdf")
        qqplot2(RS.norm.t$R, qF=function(p) sqrt(d.*qf(p, df1=d., df2=nu)),
                main.args=list(text=as.expression(substitute(bold(italic(F[list(d..,nu.)](r^2/d..))~~"Q-Q Plot"),
                    list(d..=d., nu.=nu))), side=3, cex=1.3, line=1.1, xpd=NA),
                doPDF=doPDF, file=file))

        ## Q-Q plot: R t_nu against the quantiles of F_g for a multivariate normal distribution
        file <- paste0("ggof_radial_true=t4_H0=normal_d=", d.,"_Sigma=", st, ".pdf")
        qqplot2(RS.t.norm$R, qF=function(p) sqrt(qchisq(p, df=d.)),
                main.args=list(text=as.expression(substitute(bold(italic(chi[d..])~~"Q-Q Plot"), list(d..=d.))),
                side=3, cex=1.3, line=1.1, xpd=NA),
                doPDF=doPDF, file=file))

        ## Q-Q plot: R t_nu against the correct quantiles
        file <- paste0("ggof_radial_true=t4_H0=t4_d=", d., "_Sigma=", st, ".pdf")
        qqplot2(RS.t.t$R, qF=function(p) sqrt(d.*qf(p, df1=d., df2=nu)),
                main.args=list(text=as.expression(substitute(bold(italic(F[list(d..,nu.)](r^2/d..))~~"Q-Q Plot"),
                    list(d..=d., nu.=nu))), side=3, cex=1.3, line=1.1, xpd=NA),
                doPDF=doPDF, file=file))

    }
}



### 3) Non-exchangeable models #################################################

3.1) generate data from 3rd board
3.2) do RS decomposition
3.3)

## Johanna Ziegel...
set.seed(.seed)
n <- 250 # sample size
d <- 2 # dimension
tau <- 0.5 # Kendall's tau
family <- "Gumbel" # family
cop <- getAcop(family) # define copula
R <- cop@V0(n, cop@iTau(tau)) # radial parts
E <- matrix(rexp(n*d), ncol=d) # exponentials
U <- t(E/rowSums(E))# (d,n) matrix; uniform distribution on the unit simplex
A <- matrix(c(1,-1/10,
              0.5,1), ncol=2, byrow=TRUE) # (d,d) matrix A
stopifnot(A[1,1]==1, A[2,2]==1)
X. <- A %*% U # (d,d)*(d,n)=(d,n)
X <- t(rep(R, each=d) * X.)
## pairs(pobs(X))


tau <- c(0.2, 0.6, 0.8)
r <- iTau(ellipCopula(familyTrue, df=df), tau)
P <- c(r[2], r[1], r[1], r[1], # upper triangle (without diagonal) of correlation "matrix"
             r[1], r[1], r[1],
                   r[3], r[3],
                         r[3])
cop <- ellipCopula(familyTrue, param=P, dim=d, dispstr="un", df=df)



### 4) Checking the density of R for the (not so well-working) *elliptical* case

## Problem: no non-parametric estimator for the density generator known.
## Density estimates still work comparably well.

## simulate t_4 data
n <- 1000
d <- 10
rho <- 0.5
mu <- rep(0, d)
Sigma <- outer(1:d, 1:d, FUN=function(i,j) rho^abs(i-j))
set.seed(.seed)
X.t <- rep(mu, each=n) + rmvt(n, sigma=Sigma, df=4) # multivariate t data

## (R,S) decomposition
RS.norm <- RSpobs(X.t, method="ellip", qQg=qnorm)
RS.t2   <- RSpobs(X.t, method="ellip", qQg=function(p) qt(p, df=2))
RS.t4   <- RSpobs(X.t, method="ellip", qQg=function(p) qt(p, df=4)) # true
RS.t10  <- RSpobs(X.t, method="ellip", qQg=function(p) qt(p, df=10))

## R density estimates
## (already under the 'assumed' model (due to (R,S) decomposition
##  in the *elliptical* case))
d.hat.norm <- density(RS.norm$R)
d.hat.t2   <- density(RS.t2$R)
d.hat.t4   <- density(RS.t4$R) # true
d.hat.t10  <- density(RS.t10$R)

## x range
xran <- c(4e-1, max(dens.norm$x, dens.t2$x, dens.t4$x, dens.t10$x))
x <- seq(xran[1], xran[2], length.out=513)

## corresponding theoretical R densities
d.norm <- 2*x*dchisq(x^2, df=d)
d.t2   <- 2*x*df(x^2/d, df1=d, df2=2)/d
d.t4   <- 2*x*df(x^2/d, df1=d, df2=4)/d
d.t10  <- 2*x*df(x^2/d, df1=d, df2=10)/d

## plot of density estimate based on
yran <- c(0, max(d.hat.norm$y, d.hat.t2$y, d.hat.t4$y, d.hat.t10$y,
                 d.norm, d.t2, d.t4, d.t10))
## normal
plot(dens.norm, xlim=xran, ylim=yran, lty=2, log="x",
     xlab=expression(italic(R)),
     main="R density estimates vs assumed, theoretical densities") # R density estimate
lines(x, d.norm) # 'true' if norm was correct
## t2
lines(dens.t2, col="red", lty=2) # R density estimate
lines(x, d.t2, col="red") # 'true' if t2 was correct
## t4
lines(dens.t4, col="darkgreen", lty=2) # R density estimate
lines(x, d.t4, col="darkgreen") # true
## t10
lines(dens.t10, col="blue", lty=2) # R density estimate
lines(x, d.t10, col="blue") # 'true' if t10 was correct
## legend
legend("topright", inset=0.04, bty="n", lty=rep(1:2, 4),
       col=rep(c("black", "red", "darkgreen", "blue"), each=2),
       legend=c(expression(italic(N)), expression(italic(N)),
                expression(italic(t)[2]), expression(italic(t)[2]),
                expression(italic(t)[4]~~"(true)"), expression(italic(t)[4]),
                expression(italic(t)[10]), expression(italic(t)[10])) )



### 5) Application #############################################################

##' @title -log-likelihood for t copulas
##' @param nu d.o.f. parameter
##' @param P standardized dispersion matrix
##' @param u data matrix (in [0,1]^d)
##' @return -log-likelihood for a t copula
##' @author Marius Hofert
nLLt <- function(nu, P, u) {
    stopifnot(require(mvtnorm))
    stopifnot((d <- ncol(u))==ncol(P), ncol(P)==nrow(P))
    qtu <- qt(u, df=nu)
    ldtnu <- function(u, P, nu) dmvt(qtu, sigma=P, df=nu, log=TRUE) -
        rowSums(dt(qtu, df=nu, log=TRUE)) # t copula log-density
    -sum(ldtnu(u, P=P, nu=nu))
}

## data, log-returns, pseudo-observations
data(SMI.12)
d <- ncol(x <- diff(log(SMI.12))) # log-returns
u <- pobs(x) # pseudo-observations
tau <- cor(u, method="kendall")

## TODO: visualize matrix of pairwise tau! => not Archimedean
image(M, colorkey=TRUE, col.regions=gray((1:15)/16))

## Estimate a multivariate t copula (with the approach of Demarta, McNeil (2005))
## assuming [and later checking] if the data comes from a meta-t-model

## estimate P
P <- as.matrix(nearPD(sin(cor(x, method="kendall")*pi/2), corr=TRUE)$mat)

## estimate nu via MLE for given P
nu <- seq(.5, 64, by=.25)
nLL <- sapply(nu, nLLt, P=P, u=u)
plot(nu, nLL+1200, type="l", log="xy",
     main=expression(bold("Negative log-likelihood as a function in"~~nu)),
     xlab=bquote(nu), ylab=expression(1200-logL(nu)))
## now that we got the picture, find the minimum:
nu. <- optimize(nLLt, interval=c(.5, 64), P=P, u=u, tol=1e-7)$minimum # 11.96

## Checking if the data indeed comes from a meta-t_{nu.}-model
RS.t <- RSpobs(x, method="ellip")
R.t <- RS.t$R
S.t <- RS.t$S

## TODO: new... not sure yet what we should do in this case.

## ## Q-Q plot of R against the quantiles of F_R for the estimated t distribution
## qqplot2(R.t, qF=function(p) sqrt(d*qf(p, df1=d, df2=nu.)),
##         main.args=list(text=as.expression(substitute(bold(italic(F[list(d.,nu.)](r^2/d.))~~"Q-Q Plot"),
##             list(d.=d, nu.=round(nu.,2)))), side=3, cex=1.3, line=1.1, xpd=NA))

## ## Q-Q plot of the angular distribution (Bmat[,k] should follow a Beta(k/2, (d-k)/2) distribution)
## Bmat <- gofBTstat(S.t)
## qqp(1, Bmat=Bmat) # k=1
## qqp(10, Bmat=Bmat) # k=10

## ## Check independence between radial part and B_1 and B_3
## plot(pobs(cbind(R.t, Bmat[,1])), xlab=expression(italic(R)), ylab=expression(italic(B)[1]),
##      main=expression(bold("Rank plot between"~~italic(R)~~"and"~~italic(B)[1])))
## plot(pobs(cbind(R.t, Bmat[,10])), xlab=expression(italic(R)), ylab=expression(italic(B)[10]),
##      main=expression(bold("Rank plot between"~~italic(R)~~"and"~~italic(B)[10])))
