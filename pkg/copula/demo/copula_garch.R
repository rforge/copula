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

## The copula--GARCH model

require(copula)
require(rugarch)


### 1) Simulate data from two ARMA(1,1)-GARCH(1,1) processes with dependent innovations

## simulate innovation distribution
n <- 200 # sample size
d <- 2 # dimension
family <- "t" # copula family
nu <- 3 # degrees of freedom
tau <- 0.5 # Kendall's tau
th <- iTau(ellipCopula(family, df=nu), tau) # corresponding parameter
cop <- ellipCopula(family, param=th, dim=d, df=nu) # define copula object
set.seed(271) # set seed
U <- rCopula(n, cop) # sample the copula
Z <- qnorm(U) # adjust margins

## simulate joint ARMA(1,1)-GARCH(1,1) process with these innovations
## Recall: ARMA(p_1,q_1)-GARCH(p_2,q_2) model:
##         X_t = mu_t + sigma_t * Z_t
##        mu_t = mu + \sum_{k=1}^{p_1} \phi_k*(X_{t-k}-\mu) + \sum_{k=1}^{q_1} \theta_k*(X_{t-k}-\mu_{t-k})
##   sigma_t^2 = \alpha_0 + \sum_{k=1}^{p_2} \alpha_k*(X_{t-k}-\mu_{t-k})^2 + \sum_{k=1}^{q_2} \beta_k*\sigma_{t-k}^2
## set parameters
mu <- 1
ar1 <- 0.5
ma1 <- 0.3
omega <- 2 # alpha_0 (conditional variance intercept)
alpha1 <- 0.4
beta1 <- 0.2
uspec <- ugarchspec(mean.model = list(armaOrder=c(1,1)),
                    variance.model = list(model = "sGARCH", garchOrder=c(1,1)), # standard GARCH
                    fixed.pars=list(mu=mu, ar1=ar1, ma1=ma1, omega=omega, alpha1=alpha1, beta1=beta1),
                    distribution.model = "norm") # conditional innovation density
## note: ugarchpath(): simulate from a spec; ugarchsim(): simulate from a fitted object
X <- ugarchpath(uspec, n.sim=n, # simulated path length
                m.sim=d, # number of paths to simulate
                custom.dist=list(name="sample", distfit=Z)) # pass a sample (n.sim, m.sim)-matrix
str(X, max.level=2) # => @path$sigmaSim, $seriesSim, $residSim
matplot(X@path$sigmaSim, type="l") # plot of sigma's (conditional standard deviations)
matplot(X@path$seriesSim, type="l") # plot of X's
matplot(X@path$residSim, type="l") # plot of Z's
plot(pobs(X@path$residSim)) # plot of Z's pseudo-observations => seem fine


### 2) Fit procedure based on the simulated data ###############################

## fit ARMA(1,1)-GARCH(1,1) process to X
## remove 'fixed.pars' from specification to be able to fit
uspec <- ugarchspec(mean.model = list(armaOrder=c(1,1)),
                    variance.model = list(model = "sGARCH", garchOrder=c(1,1)), # standard GARCH
                    distribution.model = "norm") # conditional innovation density
fit <- apply(X@path$seriesSim, 2, function(x) ugarchfit(uspec, x))
str(fit, max.level=3)
str(fit[[1]], max.level=2) # for first time series
stopifnot(identical(fit[[1]]@fit$residuals, residuals(fit[[1]]@fit))) # => the same

## check residuals
Z. <- sapply(fit, function(fit.) residuals(fit.@fit))
U. <- pobs(Z.)
plot(U., xlab=expression(italic(hat(U)[1])), ylab=expression(italic(hat(U)[2]))) # plot of Z's pseudo-observations => seem fine

## fit a t copula to the residuals Z
fitcop <- fitCopula(ellipCopula("t", dim=2), data=U., method="mpl")
fitcop@estimate # hat{rho}, hat{nu}; close to th, nu


### 3) Simulate from the fitted model ##########################################

## simulate from the fitted copula model
U.. <- rCopula(n, fitcop@copula)
Z.. <- qnorm(U..)

## simulate from the fitted time series model
X..sim <- lapply(1:2, function(j)
                 ugarchsim(fit[[j]], n.sim=n, m.sim=1,
                           custom.dist=list(name="sample",
                           distfit=Z..[,j, drop=FALSE]))@simulation)
str(X..sim, max.level=3)
X..Z <- sapply(X..sim, function(x) x$residSim)
plot(X..Z)
X.. <- sapply(X..sim, function(x) x$seriesSim)
matplot(pobs(X..), type="l")

