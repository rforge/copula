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


### setup ######################################################################

## load packages
stopifnot(require(grid), require(lattice), require(partitions))


### plot device ################################################################

##' Print a plot to pdf and crop it
##'
##' @title Print and crop a plot to .pdf
##' @param plot the plot, for example, a ggplot
##' @param file the file name, for example "out/Rplots.pdf"
##' @param width the width of the pdf
##' @param height the height of the pdf
##' @param crop logical indicating whether the resulting .pdf should be cropped
##' @param ... additional arguments passed on to pdf()
##' @return none
##' @author Marius Hofert
##' Note: this requires pdfcrop and pdftex being found in PATH
pdf.trellis <- function(plot, file, width=10, height=6, crop=TRUE, ...){
    trellis.device("pdf", file=file, width=width, height=height, ...) # for lattice
    print(plot)
    dev.off()
    if(crop){
        f <- file.path(getwd(), file)
        system(paste("pdfcrop --pdftexcmd pdftex", f, f, "1>/dev/null 2>&1"),
               intern=FALSE)
    }
    invisible()
}


### functions for the likelihood ###############################################

##' Compute the coefficients a_{s, d_s k}(t) for proper child copulas
##'
##' @title Compute the coefficients a_{s, d_s k}(t) for all t and k in 1:d_s
##' @param t vector t_s(u_s)
##' @param cops sector s copula list (of type nacopula with no children)
##' @param th0 parameter theta of the root copula
##' @return n x d_s matrix containing a_{s, d_s k}(t) = for all t and all k in 1:d_s
##' @author Marius Hofert
##' Note: s_{nk}(x) = \sum_{l=k}^n s(n,l)S(l,k)x^l = (-1)^{n-k} copula:::coeffG(n, x)
##'       but we need s_{nk}(x) vectorized in k and (possibly; for AMH) x
a.coeff <- function(t, cops, th0)
{
    stopifnot(is(cops, "nacopula"), length(cops@childCops)==0)
    ds <- dim(cops) # dimension d_s of the sector copula cops
    k. <- 1:ds
    ths <- cops@copula@theta # parameter theta_s of the sector copula cops
    switch(cops@copula@name, # sector copula family
           "AMH"={
               stop("a.coeff: AMH's family is currently not supported")
               ## Argh... the following code fails because copula:::coeffG does
               ## not support alpha > 1
               th0s <- (ths-th0)/(1-ths)
               n <- length(t)
               coeff <- sapply(t, function(t.) copula:::coeffG(ds,
                                                               1/(1-th0s*exp(-t.))))
               ## => coeff is an d_s x n matrix (rows = k., cols = t)
               t(coeff * (-1)^(ds-k.)) # n x d_s matrix
           },
           "Clayton"={
               as <- th0/ths
               n <- length(t)
               t. <- outer(1+t, as*k.-ds, FUN="^") # n x d_s matrix
               t. * rep((-1)^(ds-k.)*copula:::coeffG(ds, as), each=n)
           },
           "Frank"={
               stop("a.coeff: Frank's family is currently not supported")
           },
           "Gumbel"={
               as <- th0/ths
               n <- length(t)
               t. <- outer(t, as*k.-ds, FUN="^") # n x d_s matrix
               t. * rep((-1)^(ds-k.)*copula:::coeffG(ds, as), each=n)
           },
           "Joe"={
               stop("a.coeff: Joe's family is currently not supported")
           },
           stop("Family not supported"))
}

##' Coefficients b_k(t), k=d_0,..,d for all t's (sample size n) where
##'
##'     b_k(t) = \sum_{j in Q_{d.,k}^{d_0} \prod_{s=1}^{d_0} a_{s, d_s j_s}(t_s(u_s))
##'
##' where t_s is the sum over psiInv of sector s,
##'       a_{s, d_s j_s}(t) is a function specific to each family (and = 1 for
##'       non-sectorial parts),
##'   and Q_{d.,k}^{d_0} is the set of all vectors j=(j_1,..,j_d_0) such that the sum
##'       of the components j_1,..,j_d_0 equals k and j_s <= d_s for all s=1,..,d_0.
##' @title Coefficients b_k(t), k=d_0,..,d
##' @param t n x S matrix where the sth column contains the sum of psi_s^{-1}'s
##'        of the corresponding (proper) sector s
##' @param cop nacopula object of nesting depth 2 with S sectors and possibly
##'        non-sectorial parts (d_0 = dim(C0))
##' @return n x (d-d_0+1) matrix where the kth column contains b_k(t), k=d_0,..,d
##' @author Marius Hofert
b.coeff <- function(t, cop)
{
    stopifnot(t>=0, is(cop, "outer_nacopula"), nesdepth(cop)==2,
              (S <- ncol(t)) == (lchild <- length(cop@childCops)))
    if(!is.matrix(t)) t <- rbind(t) # convert t to a matrix (if necessary)

    ## get copula dimensions
    C <- cop
    d <- dim(C) # dimension of C
    d0 <- length(C@comp) + lchild # dimension of C0
    d. <- rep(1, d0) # vector of dimension d0 giving the dimensions of each sector (including the degenerate ones in the end with dimension 1 each)
    d.[1:S] <- unlist(lapply(C@childCops, dim)) # build in the dimensions of each proper child cop (in the beginning)
    stopifnot(sum(d.)==d)

    ## build the matrices Q for each k=d0,..,d (according to an email
    ## conversation with Robin Hankin = maintainer("partitions") from 2012-01-21)
    Q <- lapply(d0:d, function(k) blockparts(d.-rep(1L, d0), k-d0)+1L)
    ## => the kth element of this list Q now contains an d_0 x ? matrix where ? is
    ##    the number of elements of Q_{d.,k}^{d_0}. Each column thus contains one
    ##    vector j in Q_{d.,k}^{d_0}

    ## build the list A of length d0 containing n x d_s matrices
    A <- lapply(1:d0, function(s){
        if(d.[s]==1){ # non-sectorial components
            matrix(rep(1, n), ncol=1)
        } else { # proper child copulas
            a.coeff(t[,s], cops=C@childCops[[s]], th0=C@copula@theta)
        }
    })
    ## => the sth (in 1:d0) element of this list contains a n x d_s matrix with
    ##    the values a_{s, d_s k}, k=1,..,d_s.

    ## build b_k, k=d0,..,d
    n <- nrow(t)
    b <- matrix(, nrow=n, ncol=d-d0+1)
    for(k in d0:d){ # compute b_k=b[,k-d0+1] (first arg = t)
        ind <- k-d0+1
        Q.k <- Q[[ind]] # pick out Q for this fixed k
        ## walk over the rows of Q.k (they contain all j[s] for s fixed)
        factor.t.s.js <- array(, dim=c(n, d0, ncol(Q.k))) # contains the factors
        for(s in 1:d0){
            j.s <- Q.k[s, ] # sth row of Q.k
            factor.t.s.js[, s, ] <- A[[s]][, j.s] # n x ncol(Q.k) matrix
        }
        prods <- apply(factor.t.s.js, MARGIN=c(1,3), FUN=prod) # n x ncol(Q.k) matrix
        b[,ind] <- rowSums(prods)
    }
    b # return
}


### likelihood #################################################################

##' Log-likelihood of a two-level nested Archimedean copula
##'
##' @title Log-likelihood of a two-level nested Archimedean copula
##' @param cop two-level nested Archimedean copula ("outer_nacopula")
##' @param u matrix of realizations/observations u
##' @return -log-likelihood
##' @author Marius Hofert
nacLL <- function(cop, u)
{
    if(!is.matrix(u)) u <- rbind(u)
    stopifnot(is(cop, "outer_nacopula"), nesdepth(cop)<=2,
              (d <- ncol(u))==dim(cop),
              cop@copula@theta <= min(unlist(lapply(cop@childCops, function(cc)
              cc@copula@theta)))) # sufficient nesting condition
    ## setup
    C <- cop
    S <- length(C@childCops) # number of proper child copulas (<= d0)
    if(nesdepth(C)==1 || S == 0) stop("nacLL: You work with an Archimedean copula, so use the more efficient algorithms implemented for this special case.")
    lcomp <- length(C@comp) # if > 0 there is a non-sectorial part
    d0 <- S+lcomp # dim(C_0)
    n <- nrow(u)
    lpsiInvD1. <- matrix(, nrow=n, ncol=if(lcomp>0) S+1 else S) # n x S(+1) matrix of log(psiInvD1abs()); the "+1" comes from the non-sectorial part (if there is one)
    eta0. <- matrix(, nrow=n, ncol=d0) # n x d_0 matrix containing C_s(u_s)'s and u_0's (if there is a non-sectorial part)

    ## walk over the non-sectorial part (if available)
    if(lcomp > 0){
        u0 <- u[,C@comp, drop=FALSE]
        lpsiInvD1.[,S+1] <- rowSums(C@copula@psiInvD1abs(u0, theta=C@copula@theta,
                                                         log=TRUE))
        eta0.[,(S+1):d0] <- u0 # just set to the u's of the non-sectorial part
    }

    ## now walk over all sectors
    t. <- matrix(, nrow=n, ncol=S) # n x S matrix of t's
    for(s in 1:S){
        Cs <- C@childCops[[s]] # sector s copula list
        us <- u[, Cs@comp] # col-indices of u that belong to sector s
        t.[,s] <- rowSums(Cs@copula@psiInv(us, theta=Cs@copula@theta))
        lpsiInvD1.[,s] <- rowSums(Cs@copula@psiInvD1abs(us, theta=Cs@copula@theta,
                                                        log=TRUE))
        eta0.[,s] <- Cs@copula@psi(t.[,s], theta=Cs@copula@theta) # C_s(u_s)
    }

    ## finish computations
    lpsiInvD1sum <- rowSums(lpsiInvD1.) # sum(log(psiInvD1abs)); vector of length n
    eta0 <- rowSums(C@copula@psiInv(eta0., theta=C@copula@theta)) # eta0; vector of length n; = C@copula@psiInv(pnacopula(C, u), theta=C@copula@theta)

    ## compute b_k's, k=d0,..,d
    ## note: it suffices to give b.coeff only the sectorial t's
    b.mat <- b.coeff(t., cop=C) # n x (d-d_0+1) matrix

    ## compute x_k, k=d0,..,d, a n x (d-d_0+1) matrix
    th0 <- C@copula@theta
    x <- sapply(d0:d, function(k){
        log((-1)^(d-k) * b.mat[,k-d0+1]) +
            C@copula@psiDabs(eta0, theta=th0, degree=k, log=TRUE)
    })

    ## return
    sum(copula:::lsum(t(x)) + lpsiInvD1sum) # sum over all t's
}


### Example 1: ((1,2), (3,4,5))-Gumbel #########################################

## setup
n <- 250
family <- "Gumbel"
cop. <- getAcop(family)
tau <- c(0.2, 0.4, 0.6)
th <- cop.@tauInv(tau)
cop <- onacopulaL(family, list(th[1], NULL, list(list(th[2], 1:2),
                                                 list(th[3], 3:5))))

## sample and compute log-likelihood
U <- rnacopula(n, cop)
nacLL(cop, u=U) # log-likelihood at correct parameters


### Example 2: (1, (2,3), 4, (5,6,7)) Gumbel ###################################

## setup
n <- 250
family <- "Gumbel"
cop. <- getAcop(family)
tau <- c(0.2, 0.4, 0.6)
th <- cop.@tauInv(tau)
cop <- onacopulaL(family, list(th[1], c(1,4), list(list(th[2], 2:3), list(th[3],
                                                                          5:7))))

## Sample and compute log-likelihood
U <- rnacopula(n, cop)
nacLL(cop, u=U) # log-likelihood at correct parameters


### Example 3: (1, (2,3)) Gumbel ###############################################

## setup
n <- 250
family <- "Gumbel"
cop. <- getAcop(family)
tau <- c(0.25, 0.5)
th <- cop.@tauInv(tau)
cop <- onacopulaL(family, list(th[1], 1, list(list(th[2], 2:3))))

## Sample and compute log-likelihood
U <- rnacopula(n, cop)
nacLL(cop, u=U) # log-likelihood at correct parameters


### -log-Likelihood plots ######################################################

### (1) (1, (2,...)) structure #################################################

## setup: use this to adjust the family and scomp
comp <- 1 # non-sectorial indices
stopifnot(comp==1) # otherwise, sub (see below) is wrong
scomp <- 2:3 # 2:10 or 2:3 (sectorial indices)
family <- "Gumbel" # "Clayton" or "Gumbel"

## determine values
n <- 100
cop. <- getAcop(family)
tau <- c(0.25, 0.5)
th <- cop.@tauInv(tau)
h <- 0.2
th0 <- cop.@tauInv(c(tau[1]-h, tau[1]+h))
th1 <- cop.@tauInv(c(tau[2]-h, tau[2]+h))
m <- 20 # number of grid points
th0. <- seq(th0[1], th0[2], length.out=m)
th1. <- seq(th1[1], th1[2], length.out=m)
grid <- expand.grid(th0=th0., th1=th1.)
cop <- onacopulaL(family, list(th[1], comp, list(list(th[2], scomp)))) # copula
U <- rnacopula(n, cop) # sample
nLL <- function(th, u, family, comp, scomp){
    if(th[1] > th[2]) return(NA) # sufficient nesting condition
    cop <- onacopulaL(family, list(th[1], comp, list(list(th[2], scomp))))
    -nacLL(cop, u=u)
}
val.grid <- apply(grid, 1, function(x) nLL(x, u=U, family=family, comp=comp,
                                           scomp=scomp))

## determine plot supplements
true.val <- c(th0=th[1], th1=th[2],
              nLL=nLL(th, u=U, family=family, comp=comp, scomp=scomp)) # true value
ind <- which.min(val.grid)
opt.val <- c(grid[ind,], nLL=val.grid[ind]) # optimum on the grid
pts <- rbind(true.val, opt.val) # points to add to wireframe plot
title <- paste("-log-likelihood of a nested", family, "copula") # title
mysec <- if(length(scomp)==2) bquote(italic(u[3])) else
substitute(list(...,italic(u[j])), list(j=max(scomp)))
sub <- substitute(italic(C(bolditalic(u)))==italic(C[0](u[1],C[1](u[2],mysec.)))
                  ~~~~~~ italic(n)==n. ~~~~~~ tau(theta[0])==tau0. ~~~~~~
                  tau(theta[1])==tau1., list(mysec.=mysec, n.=n, tau0.=tau[1],
                     tau1.=tau[2]))
sub <- as.expression(sub) # lattice "bug" (only needed by lattice)
xlab <- expression(italic(theta[0]))
ylab <- expression(italic(theta[1]))
zlab <- list(as.expression(-log~L*group("(",italic(theta[0])*"," ~
    italic(theta[1])*";"~bolditalic(u),")")), rot=90)

## wireframe plot
(wf <- wireframe(val.grid~grid[,1]*grid[,2], aspect=1, zoom=1.02, xlim=th0, ylim=th1,
                 zlim=c(min(val.grid, as.numeric(pts[,3]), na.rm=TRUE),
                 max(val.grid, as.numeric(pts[,3]), na.rm=TRUE)),
                 xlab=xlab, ylab=ylab, zlab=zlab, main=title, sub=sub, pts=pts,
                 par.settings=list(standard.theme(color=FALSE),
                 layout.heights=list(sub=2.4), background=list(col="#ffffff00"),
                 axis.line=list(col="transparent"), clip=list(panel="off")),
                 alpha.regions=0.5, scales=list(col=1, arrows=FALSE),
                 ## add wire/points
                 panel.3d.wireframe=function(x, y, z, xlim, ylim, zlim, xlim.scaled,
                 ylim.scaled, zlim.scaled, pts, ...){
                     panel.3dwire(x=x, y=y, z=z, xlim=xlim, ylim=ylim, zlim=zlim,
                                  xlim.scaled=xlim.scaled, ylim.scaled=ylim.scaled,
                                  zlim.scaled=zlim.scaled, ...)
                     panel.3dscatter(x=as.numeric(pts[,1]), y=as.numeric(pts[,2]),
                                     z=as.numeric(pts[,3]), xlim=xlim, ylim=ylim, zlim=zlim,
                                     xlim.scaled=xlim.scaled, ylim.scaled=ylim.scaled,
                                     zlim.scaled=zlim.scaled, type="p", col=1,
                                     pch=c(3,4), lex=2, cex=1.4, .scale=TRUE, ...)
                 }, key=list(x=-0.01, y=1, points=list(pch=c(3,4), col=1, lwd=2,
                                           cex=1.4),
                    text=list(c(expression(group("(",list(theta[0],theta[1]),")")^T),
                    expression(group("(",list(hat(theta)["0,n"],hat(theta)["1,n"]),")")^T))),
                    padding.text=3, cex=1, align=TRUE, transparent=TRUE)) )
## pdf.trellis(wf, file=paste("wf_nLL_",family,"_d=",dim(cop),".pdf",sep=""),
##             width=6, height=6)

## levelplot
xlim. <- c(min(grid[,1]), max(grid[,1]))
ylim. <- c(min(grid[,2]), max(grid[,2]))
xeps <- (xlim.[2] - xlim.[1]) * 0.04
yeps <- (ylim.[2] - ylim.[1]) * 0.04
(lp <- levelplot(val.grid~grid[,1]*grid[,2], aspect=1, xlab=xlab, ylab=ylab,
                 par.settings=list(layout.heights=list(main=3, sub=2),
                 regions=list(col=gray(140:400/400))),
                 xlim=c(xlim.[1]-xeps, xlim.[2]+xeps),
                 ylim=c(ylim.[1]-yeps, ylim.[2]+yeps),
                 main=title, sub=sub, pts=pts,
                 scales=list(alternating=c(1,1), tck=c(1,0)), contour=TRUE,
                 panel=function(x, y, z, pts, ...){
                     panel.levelplot(x=x, y=y, z=z, ...)
                     grid.points(x=pts[1,1], y=pts[1,2], pch=3,
                                 gp=gpar(lwd=2, col="black")) # + true value
                     grid.points(x=pts[2,1], y=pts[2,2], pch=4,
                                 gp=gpar(lwd=2, col="black")) # x optimum
                 }, key=list(x=0.18, y=1.09, points=list(pch=c(3,4), col=1,
                                             lwd=2, cex=1.4), columns=2,
                    text=list(c(expression(group("(",list(theta[0],theta[1]),")")^T),
                    expression(group("(",list(hat(theta)["0,n"],hat(theta)["1,n"]),")")^T))),
                    align=TRUE, transparent=TRUE)) )
## pdf.trellis(lp, file=paste("lp_nLL_",family,"_d=",dim(cop),".pdf",sep=""),
##             width=6, height=6)

