## Copyright (C) 2010 Marius Hofert and Martin Maechler
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

library(nacopula) 
options(warn = 1) 

#### ==== Demo for working with the two-parameter Archimedean GIG copula =======

source(system.file("Rsource", "GIG.R", package="nacopula"))

## ==== setup ==================================================================

## set seed
set.seed(1) 

## ==== specify parameters ====

d <- 5 # dimension
tau <- 0.2 # specify Kendall's tau
## specify the theta vector such that Kendall's tau equals the given tau
theta <- tauInv.GIG(tau, theta=c(1, NA))

## initial interval bounds for optimization
h <- c(0.3, 0.3) # for initial interval for theta_1 and theta_2
I <- matrix(c(theta[1]-h[1], theta[1]+h[1], theta[2]-h[2], theta[2]+h[2]), ncol=2,
            dimnames=list(c("lower", "upper"), c("nu", "theta"))) # initial intervals for theta_1 and theta_2
rb <- min((I[2,1]-I[1,1])/2, (I[2,2]-I[1,2])/2) # for rhobeg

## ==== sample the copula ====

n <- 1000
U <- rnacopula.GIG(n, d, theta)
tau.hat <- cor(U, method="kendall")
tau.hats <- tau.hat[upper.tri(tau.hat)]

## plot the sample
pairs(U, gap=0)

## ==== estimate the parameters ====

n <- 100 # choose small sample size due to run time [depending on psi^{-1}]
U <- rnacopula.GIG(n, d, theta)

## choose random initial point in a reasonable interval
start <- c(runif(1, min=I[1,1], max=I[2,1]), runif(1, min=I[1,2], max=I[2,2]))

## MLE
system.time(res <- optimx(par=start, method="bobyqa",
                          fn=function(x) -sum(dacopula.GIG(U, theta=x, n.MC=0, log=TRUE)), 
                          lower=c(I[1,1],I[1,2]), upper=c(I[2,1],I[2,2]), hessian=FALSE,
                          control=list(rhobeg=rb)))
res
## note: run time is mainly determined by the evaluations of psiInv.GIG

## ==== plot Kendall's tau as a function in theta for different nu's ====

th1 <- c(0, 0.003, 0.5, 1, 5, 10)
cols <- colorRampPalette(c("red", "orange", "darkgreen", "turquoise", "blue"), 
                         space="Lab")(length(th1))
for(i in seq_along(th1))
    curve(tau.GIG(cbind(th1[i],x)), 1e-12, 2, 
          main="Kendall's tau for the GIG family", ylim=c(0,1),
          xlab=expression(theta), ylab=expression(tau(nu,theta)), add=(i>1), 
          lwd=1.4, col=cols[i])
label <- as.expression(lapply(1:length(th1), function(i) substitute(nu==nu., list(nu.=th1[i]))))
legend("topright", label, bty="n", lwd=1.4, col=cols)
