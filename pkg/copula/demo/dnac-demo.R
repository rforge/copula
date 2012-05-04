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


### Densities of two-level nested Archimedean copulas ##########################

## Note: see http://arxiv.org/abs/1204.2410 for more details


### Setup ######################################################################

source(system.file("Rsource", "dnac.R", package="copula"))

require(copula)
require(grid)
require(lattice)


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
