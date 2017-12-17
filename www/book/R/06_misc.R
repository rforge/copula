## By Marius Hofert, Ivan Kojadinovic, Martin Maechler, Jun Yan

## R script for Chapter 6 of Elements of Copula Modeling with R


library(copula)


### 6.1 Ties ###################################################################

### Computing ranks in the presence of ties

set.seed(1979)
(U <- runif(8))

R.avg    <- rank(U) # ties.method = "average"
R.random <- rank(U, ties.method = "random")
R.max    <- rank(U, ties.method = "max")
R.min    <- rank(U, ties.method = "min")
stopifnot(R.random == R.avg, R.max == R.avg, R.min == R.avg)

b <- 10 # number of bins
(U.ties <- cut(U, breaks = 0:b/b, labels = 0.5:(b - 0.5)/b)) # a factor

rank(U.ties) # ties.method = "average"

rank(U.ties, ties.method = "max") # maximum rank for each tied observation

rank(U.ties, ties.method = "min") # minimum rank for each tied observation

set.seed(8)
rank(U.ties, ties.method = "random")

set.seed(46)
rank(U.ties, ties.method = "random")


### Effect of ties on multivariate scaled ranks

n <- 200
set.seed(4809)
U <- rCopula(n, claytonCopula(5))
b <- 10 # number of bins
U.ties <- (apply(U, 2, cut, breaks = 0:b/b, labels = FALSE) - 0.5) / b

ties.method <- "max" # can be changed to "min", "average" or "random"
stopifnot(all.equal(pobs(U.ties, ties.method = ties.method),
                    apply(U.ties, 2, rank, ties.method = ties.method) /
                    (nrow(U.ties) + 1))) # check

plot(pobs(U), xlab = "", ylab = "")
plot(pobs(U.ties), xlim = 0:1, ylim = 0:1, xlab = "", ylab = "")
set.seed(732)
plot(pobs(U.ties, ties.method = "random"), xlab = "", ylab = "")

a <- runif(1, min = 0, max = 1/(2 * b)) # or a fixed number in (0, 1/(2b))
set.seed(732)
V <- pobs(U.ties + runif(2 * n, min = -a, max = a))
set.seed(732)
stopifnot(all.equal(V, pobs(U.ties, ties.method = "random"))) # check


### Estimation of copula parameters in the presence of ties

theta <- iTau(frankCopula(), tau = 0.25) # copula parameter
fc25 <- frankCopula(theta) # corresponding copula object

## Discretizes all components of U; returns a matrix of real numbers
discrAll <- function(U, b)
    (apply(U, 2L, cut, breaks = 0:b/b, labels = FALSE) - 1/2) / b

##' @title Breaking the ties m times and averaging the estimates
##' @param m number of replications
##' @param U.ties discretized sample
##' @param cop copula to be fitted
##' @param method fitting method of fitCopula()
##' @param optim.method optimization method of fitCopula()
##' @return average estimate over randomly broken ties when computing pobs()
fitCopulaRand <- function(m, U.ties, cop, method, optim.method)
{
    fit1 <- function() {
        V <- pobs(U.ties, ties.method = "random") # break ties at random
        fitCopula(cop, data = V, method = method,
                  optim.method = optim.method)@estimate # return param. est.
    }
    mean(replicate(m, fit1())) # average of estimates over m randomizations
}

##' @title Parameter estimation from discrete data with various methods
##' @param n sample size
##' @param cop one-parameter copula object
##' @param discretize discretizing function
##' @param b number of bins
##' @param m number of randomizations
##' @param optim.method MPLE optimization method
##' @return parameter estimates for the different methods
oneFit <- function(n, cop, discretize, b = 10, m = 30,
                   optim.method = "BFGS")
{
    U <- rCopula(n, copula = cop) # a sample of size n from cop
    ## Rank-based fitting from the continuous sample
    V <- pobs(U) # pseudo-observations (no ties)
    mpl <- fitCopula(cop, data = V,
                     optim.method = optim.method)@estimate # MPLE
    itau <- fitCopula(cop, data = V,
                      method = "itau")@estimate # tau inversion
    ## Fitting from the discretized sample based on average ranks
    U.ties <- discretize(U, b = b) # discretize the sample
    W <- pobs(U.ties) # corresponding multivariate scaled average ranks
    mpl.ave  <- fitCopula(cop, data = W,
                          optim.method = optim.method)@estimate # MPLE
    itau.ave <- fitCopula(cop, data = W,
                          method = "itau")@estimate # tau inversion
    ## Average fits over randomized ranks based on the discretized sample
    mpl.rand  <- fitCopulaRand(m, U.ties, cop = cop, method = "mpl",
                               optim.method = optim.method) # MPLE
    itau.rand <- fitCopulaRand(m, U.ties, cop = cop, method = "itau",
                               optim.method = optim.method) # tau inversion
    ## Return the different estimates of theta
    c(mpl = mpl, mpl.ave = mpl.ave, mpl.rand = mpl.rand,
      itau = itau, itau.ave = itau.ave, itau.rand = itau.rand)
}

set.seed(2010)
est.fc25 <- withTime(replicate(100, oneFit(n = 100, cop = fc25,
                                           discretize = discrAll)))

##' @title Relative bias, standard deviation and root mean squared error
##' @param est simulation object
##' @param cop underlying copula object
##' @return rel. bias, standard dev. and root mean squared error (in %)
sumSims <- function(est, cop)
{
    tau <- tau(cop) # population version of Kendall's tau
    cat(describeCop(cop, kind = "very short"),
        "with a Kendall's tau of", tau, "\n")
    print(est$sys.time) # print user time (first component of system.time())
    ## Estimates on Kendall's tau scale
    tau.n.s <- apply(est$value, 1:2, function(x) tau(setTheta(cop, x)))
    ## Relative biases on Kendall's tau scale
    bias <- (rowMeans(tau.n.s) - tau) / tau
    ## Standard deviations of estimates on Kendall's tau scale
    std <- apply(tau.n.s, 1, sd)
    ## Root mean squared errors on Kendall's tau scale
    rmse <- sqrt(rowMeans((tau.n.s - tau)^2))
    round(rbind(bias = bias, std = std, rmse = rmse) * 100, 2)
}
sumSims(est = est.fc25, cop = fc25)

fc5 <- frankCopula(iTau(frankCopula(), tau = 0.5))
set.seed(2010)
est.fc5 <- withTime(replicate(100, oneFit(n = 100, cop = fc5,
                                          discretize = discrAll)))
sumSims(est = est.fc5, cop = fc5)

fc75 <- frankCopula(iTau(frankCopula(), tau = 0.75))
set.seed(2010)
est.fc75 <- withTime(replicate(100, oneFit(n = 100, cop = fc75,
                                           discretize = discrAll)))
sumSims(est = est.fc75, cop = fc75)

## An alternative definition of the discretizing function
## producing ties only in the first component sample
discrFirst <- function(U, b)
    cbind(cut(U[,1], breaks = 0:b/b, labels = 0.5:(b - 0.5)/b), U[,2])


### Tests not adapted for ties

##' @title Auxiliary function for computing the empirical levels of tests
##'        of exchangeability and extreme-value dependence, and
##'        parametric bootstrap and multiplier goodness-of-fit tests
##' @param n sample size
##' @param cop copula object (from which the data is generated)
##' @param discretize discretization function
##' @param b number of bins
##' @return p-values (numeric(4))
pvalTies <- function(n, cop, discretize, b)
{
    U.ties <- discretize(rCopula(n, cop), b = b) # binned samples
    c(exch = exchTest(U.ties, ties = FALSE)$p.value,
      ev   = evTestC(U.ties)$p.value,
      pb   = gofCopula(cop, U.ties, optim.method = "BFGS",
                          ties = FALSE)$p.value,
      mult = gofCopula(cop, U.ties, optim.method = "BFGS",
                          sim = "mult")$p.value)
}

gc <- gumbelCopula(iTau(gumbelCopula(), tau = 0.5)) # Gumbel-Hougaard copula
set.seed(4478)
pv <- withTime(replicate(100, pvalTies(n = 100, cop = gc,
                                       discretize = discrAll, b = 10)))

pv$sys.time # user time

alpha <- c(0.01, 0.05, 0.1) # nominal levels
rbind(nom.level = alpha,
      emp.level.exch = ecdf(pv$value["exch",])(alpha),
      emp.level.ev   = ecdf(pv$value["ev",  ])(alpha),
      emp.level.pb   = ecdf(pv$value["pb",  ])(alpha),
      emp.level.mult = ecdf(pv$value["mult",])(alpha))

summary(pv$value["mult",])


### Parametric bootstrap-based goodness-of-fit test adapted for ties

##' @title Auxiliary function for computing the empirical levels of the
##'        parametric bootstrap-based goodness-of-fit test
##' @param n sample size
##' @param cop copula object (from which the data is generated)
##' @param discretize discretization function
##' @param b number of bins
##' @return p-value (numeric(1))
pvalPB <- function(n, cop, discretize, b)
{
    U <- rCopula(n, copula = cop)
    U.ties <- discretize(U, b = b)
    gofCopula(cop, x = U.ties, optim.method = "BFGS", ties = TRUE)$p.value
}

set.seed(8848)
pvAll <- withTime(replicate(100, pvalPB(n = 100, cop = gc,
                                        discretize = discrAll, b = 10)))
pvAll$sys.time # run time

set.seed(3298)
pvFirst <- withTime(replicate(100, pvalPB(n = 100, cop = gc,
                                          discretize = discrFirst, b = 10)))
pvFirst$sys.time # run time

alpha <- c(0.01, 0.05, 0.1) # nominal levels
rbind(nom.level = alpha, emp.level.all = ecdf(pvAll$value)(alpha),
      emp.level.first = ecdf(pvFirst$value)(alpha))


### Effect of ties on cross-validation

##' @title 10-fold cross-validation in the presence of ties
##' @param n sample size
##' @param cop copula object (from which the data is generated)
##' @param discretize discretization function
##' @param b number of bins
##' @param optim.method optimization method for MPLE
##' @return label of the copula with the highest xv score
xv <- function(n, cop, discretize = discrAll, b, optim.method = "BFGS")
{
    U <- rCopula(n, copula = cop)
    U.ties <- discretize(U, b = b)
    score <- c(xvCopula(claytonCopula(), x = U.ties, k = 10,
                        optim.method = optim.method),
               xvCopula(gumbelCopula(), x = U.ties, k = 10,
                        optim.method = optim.method),
               xvCopula(frankCopula(), x = U.ties, k = 10,
                        optim.method = optim.method))
    c("C", "GH", "F")[which(score == max(score))]
}

cc <- claytonCopula(iTau(claytonCopula(), tau = 0.5))
set.seed(2885)
withTime(table(replicate(100, xv(n = 100, cop = cc, b = 20))))

set.seed(2885)
withTime(table(replicate(100, xv(n = 100, cop = gc, b = 20))))

fc <- frankCopula(iTau(frankCopula(), tau = 0.5))
set.seed(2885)
withTime(table(replicate(100, xv(n = 100, cop = fc, b = 20))))


### Analysis of the loss insurance data

data(loss)
X <- as.matrix(subset(loss, censored == 0, select = c("loss", "alae")))

## Percentages of ties
100 * apply(X, 2, function(x) 1 - length(unique(x))/length(x))

U <- pobs(X)
plot(U, xlab = quote(U[1]~~"(Loss)"), ylab = quote(U[2]~~"(ALAE)"))
Y <- qnorm(U)
plot(Y, xlab = quote(Y[1]~~"(Loss)"), ylab = quote(Y[2]~~"(ALAE)"))

set.seed(3070)
withTime(round(c(  exchTest(X, ties = TRUE)$p.value,
                 radSymTest(X, ties = TRUE)$p.value,
                    evTestK(X, ties = TRUE)$p.value), 4))

## Goodness-of-fit testing
set.seed(4634)
optim.method <- "BFGS" # the numerical optimization method for MPLE
withTime(
    round(c(gofCopula(gumbelCopula(), x = X, ties = TRUE,
                      optim.method = optim.method)$p.value,
            gofCopula(rotCopula(claytonCopula()), x = X, ties = TRUE,
                      optim.method = optim.method)$p.value,
            gofCopula(frankCopula(), x = X, ties = TRUE,
                      optim.method = optim.method)$p.value,
            gofCopula(plackettCopula(), x = X, ties = TRUE,
                      optim.method = optim.method)$p.value,
            gofCopula(normalCopula(), x = X, ties = TRUE,
                      optim.method = optim.method)$p.value), 4)
)

## Model selection
set.seed(4807)
k <- 50 # for k-fold cross-validation
withTime(
    round(c(xvCopula(gumbelCopula(), x = X, k = k,
                     optim.method = optim.method),
            xvCopula(rotCopula(claytonCopula()), x = X, k = k,
                     optim.method = optim.method),
            xvCopula(frankCopula(), x = X, k = k,
                     optim.method = optim.method),
            xvCopula(plackettCopula(), x = X, k = k,
                     optim.method = optim.method),
            xvCopula(normalCopula(), x = X, k = k,
                     optim.method = optim.method)), 1)
)

withTime(exchTest(X, ties = FALSE))


### 6.2 Selected copula tests and models for time series #######################

### 6.2.1 Tests of stationarity ################################################

### Test of stationarity based on S_n^H

data(rdj)
library(xts)
Xrdj <- xts(rdj[,-1], order.by = rdj[,1])

plot.zoo(Xrdj, main = "", xlab = "", mar = c(0, 7, 0, 2.1))

library(npcp)
set.seed(981)
(res <- withTime(cpDist(Xrdj, b = NULL)))

out <- res$value # the output of the test (object of class 'htest')
rdj[which(out$cvm == out$statistic), 1]

plot(out$cvm, type = "l", xlab = "k", ylab = quote({S^H}[list(n,k)]))

data(gasoil) # oil and gas prices
library(qrmtools) # for returns()
Rgasoil <- returns(gasoil[,-1]) # bivariate daily log-returns
Xgasoil <- xts(Rgasoil, order.by = gasoil[-1, 1]) # corresponding xts object

plot.zoo(Xgasoil, main = "", xlab = "", mar = c(0, 7, 0, 2.1))

set.seed(292)
withTime(cpDist(Xgasoil, b = NULL))


### Test of stationarity based on S_n^C

set.seed(314)
withTime(cpCopula(Xrdj, b = NULL, method = "nonseq"))

set.seed(137)
withTime(cpCopula(Xgasoil, b = NULL, method = "nonseq"))


### Test of stationarity based on S_n^{C^s}

set.seed(3355)
withTime(cpAutocop(Xrdj[,1], lag = 1))

set.seed(3355)
withTime(cpAutocop(Xrdj[,1], lag = 2))

set.seed(3355)
withTime(cpAutocop(Xrdj[,1], lag = 3))

set.seed(3105)
withTime(cpAutocop(Xgasoil[,1], lag = 1))

set.seed(2895)
withTime(cpAutocop(Xgasoil[,2], lag = 1))


### 6.2.2 Tests of serial independence #########################################

### Correlogram and Ljung-Box test of serial independence

colnames(Xgasoil) <- c("Oil", "Gas")
acf(Xgasoil^2, ci.col = 1)

Box.test(Xgasoil[,1]^2, lag =  5, type = "Ljung-Box")

Box.test(Xgasoil[,2]^2, lag =  5, type = "Ljung-Box")

Box.test(Xgasoil[,1]^2, lag = 20, type = "Ljung-Box")

Box.test(Xgasoil[,2]^2, lag = 20, type = "Ljung-Box")


### Ljung-Box tests can be too liberal

pvalBox <- function(n)
{
    x2 <- rt(n, df = 4)^2
    c(lag5  = Box.test(x2, type = "Ljung-Box", lag =  5)$p.value,
      lag20 = Box.test(x2, type = "Ljung-Box", lag = 20)$p.value)
}

set.seed(3298)
pv <- replicate(10000, pvalBox(500))
alpha <- c(0.01, 0.05, 0.1) # nominal levels
rbind(nom.level = alpha,
      emp.level.lag5  = ecdf(pv["lag5", ])(alpha),
      emp.level.lag20 = ecdf(pv["lag20",])(alpha))


### Tests of serial independence based on S_n^{Pi^s}

set.seed(137)
sI.d <- withTime(serialIndepTestSim(nrow(Xgasoil), lag.max = 5))
sI.d$sys.time # the run time

serialIndepTest(Xgasoil[,1]^2, d = sI.d$value)

serialIndepTest(Xgasoil[,2]^2, d = sI.d$value)


### 6.2.3 Models for multivariate time series based on conditional copulas #####

### Conditional modeling based on ARMA--GARCH marginal models

library(rugarch)
## Specify ARMA(1,1)-GARCH(1,1) model with Student t innovations
meanModel <- list(armaOrder = c(1,1)) # ARMA(1,1)
varModel  <- list(model = "sGARCH", garchOrder = c(1,1)) # GARCH(1,1)
uspec <- ugarchspec(varModel, mean.model = meanModel,
                    distribution.model = "std") # scaled Student t

## Fit marginal ARMA-GARCH models
fit <- apply(Xgasoil, 2, function(x) ugarchfit(uspec, data = x))
## Extract the estimated standardized residuals
eps <- sapply(fit, residuals, standardize = TRUE) # standardized residuals
(nus <- sapply(fit, function(x) x@fit$coef[["shape"]])) # fitted d.o.f.

set.seed(2013)
withTime(cpCopula(eps, b = 1, method = "nonseq"))

U <- pobs(eps) # pseudo-observations from the residuals eps
plot(U, xlab = expression(U[1]), ylab = expression(U[2]))

## Fit a t copula to the estimated standardized residuals
fitcop <- fitCopula(tCopula(), data = U, method = "mpl")
fitcop@estimate # estimated correlation parameter rho and d.o.f. nu

cop <- fitcop@copula # fitted t copula

## Simulate from the bivariate model
## 1) Simulate from the fitted copula
set.seed(271) # set seed
n.sim <- 260 # sample size
m.sim <- 1000 # number of paths
U. <- rCopula(n.sim * m.sim, cop) # simulate from the fitted copula
## 2) Quantile-transform the corresponding innovations
##    Note: eps have to be standardized (mean 0, variance 1) for ugarchsim()
eps. <- sapply(1:2, function(j)
    sqrt((nus[j]-2)/nus[j]) * qt(U.[,j], df = nus[j]))
## 3) Feed the (cross-sec. dependent) innovations to the marginal ARMA-GARCH
##    models and simulate from them
sim <- lapply(1:2, function(j)
    ugarchsim(fit[[j]], # fitted marginal ARMA-GARCH model
              n.sim = n.sim, # sample size
              m.sim = m.sim, # number of trajectories/paths
              custom.dist = list(name = "sample", # our innovations
                                 distfit = matrix(eps.[,j], ncol = m.sim))))
## 4) Extract the simulated (cross-sec. dependent) series X_t and build the
##    corresponding (predicted/simulated) oil and gas prices
X. <- lapply(sim, function(x) fitted(x)) # equal to seriesSim in x@simulation
S.t <- as.numeric(tail(gasoil, n = 1)[,2:3]) # last available prices S_t
library(qrmtools) # for returns()
S. <- lapply(1:2, function(j) # predicted prices for each stock
    returns(X.[[j]], inverse = TRUE, start = rep(S.t[j], m.sim)))
S.T <- sapply(1:2, function(j) tail(S.[[j]], n = 1)) # pick out prices at T

library(MASS)
pred.dens <- kde2d(S.T[,1], S.T[,2], n = 300, lims = c(0, 200, 0, 15))
image(pred.dens, xlab = "Oil price", ylab = "Gas price",
      col = gray(seq(1, 0, length.out = 100)))


### 6.3 Regression #############################################################

### Conditional modeling based on marginal gamma GLMs

data(NELS88, package = "copulaData")
nels <- subset(NELS88, select = -ID) # remove school ID

nels$Size <- scale(nels$Size) # mean 0, standard deviation 1

##' @title Marginal conditional negative log-likelihood
##' @param beta.m parameter vector defining the marginal calibration map
##' @param x vector of values of one of the three scores
##' @param z design matrix
##' @param pobs logical indicating whether, additionally, the parametric
##'        pseudo-observations shall be computed and returned
##' @return -log-likelihood and, possibly, the parametric pseudo-observations
nmLL <- function(beta.m, x, z, pobs = FALSE)
{
    p <- ncol(z) + 1 # number of parameters
    mu.z <- exp(z %*% beta.m[1:(p-1)])
    a.z <- 1 / beta.m[p] # shape
    s.z <- mu.z * beta.m[p] # scale
    nLL <- -sum(dgamma(x, shape = a.z, scale = s.z, log = TRUE))
    if (!pobs) nLL else
        list(nLL = nLL, U = pgamma(x, shape = a.z, scale = s.z))
}

## Build the design matrix
z <- model.matrix(~ Minority + SES + Female + Public + Size +
                      Urban + Rural, data = nels)
p <- ncol(z) + 1 # number of parameters per margin

math.glm <- glm(Math ~ Minority + SES + Female + Public + Size +
                Urban + Rural, data = nels, family = Gamma(link = "log"))
(math.summary <- summary(math.glm))

## The estimates of (beta_{1,1}, ..., beta_{1,9})
(ts.math <- c(math.glm$coefficients, disp = math.summary$dispersion))

## Minimizing the marginal conditional negative log-likelihood
## using the previously obtained estimates as initial values
res <- optim(ts.math, nmLL, x = nels[,"Math"], z = z, method = "BFGS")
## Compare GLM and ML estimates: change is small
stopifnot(all.equal(ts.math, res$par, tolerance = 1e-3))

## Science score
sci.glm <- glm(Science ~ Minority + SES + Female + Public + Size +
               Urban + Rural, data = nels, family = Gamma(link = "log"))
## The estimates of (beta_{2,1}, ..., beta_{2,9})
(ts.sci <- c(sci.glm$coefficients, disp = summary(sci.glm)$dispersion))

## Reading score
read.glm <- glm(Reading ~ Minority + SES + Female + Public + Size +
                Urban + Rural, data = nels, family = Gamma(link = "log"))
## The estimates of (beta_{3,1}, ..., beta_{3,9})
(ts.read <- c(read.glm$coefficients, disp = summary(read.glm)$dispersion))

## Parametric pseudo-observations from the underlying trivariate conditional
## copula under the parametric assumptions and the simplifying assumption
V <- cbind("V[1]" = nmLL(ts.math, nels[,"Math"],    z, pobs = TRUE)$U,
           "V[2]" = nmLL(ts.sci,  nels[,"Science"], z, pobs = TRUE)$U,
           "V[3]" = nmLL(ts.read, nels[,"Reading"], z, pobs = TRUE)$U)

splom2(V, cex = 0.3, col.mat = "black") # scatter-plot matrix
library(lattice)
cloud2(V) # 3d cloud plot based on lattice's cloud()

U <- pobs(cbind(residuals(math.glm), residuals(sci.glm),
                residuals(read.glm)))

stopifnot(all.equal(U, pobs(V), check.attributes = FALSE))

summary(ts.g <- fitCopula(gumbelCopula(dim = 3), data = U))

summary(ts.f <- fitCopula(frankCopula (dim = 3), data = U))

summary(ts.n.ex <- fitCopula(normalCopula(dim = 3, dispstr = "ex"), data=U))

summary(ts.n.un <- fitCopula(normalCopula(dim = 3, dispstr = "un"), data=U))

k <- 100 # for k-fold cross-validation
optim.method <- "BFGS" # the numerical optimization method for MPLE
set.seed(3090)
withTime(
    round(c(xvCopula(gumbelCopula(dim = 3),
                     x = V, k = k, optim.method = optim.method),
            xvCopula(frankCopula(dim = 3),
                     x = V, k = k, optim.method = optim.method),
            xvCopula(normalCopula(dim = 3), # homogeneous (exchangeable)
                     x = V, k = k, optim.method = optim.method),
            xvCopula(normalCopula(dim = 3, dispstr = "un"), # unstruct. cor.
                     x = V, k = k, optim.method = optim.method)), 1)
)

set.seed(2723)
withTime(gofCopula(normalCopula(dim = 3), x = V, simulation = "mult"))

##' @title Full conditional negative log-likelihood function
##' @param par param. vector defining the marg. and copula calibration maps
##' @param x matrix of values of the responses
##' @param z design matrix
##' @param copula a trivariate one-parameter copula object
##' @return -log-likelihood
nfLL <- function(par, x, z, copula)
{
    beta <- par[1] # copula parameter
    tc <- tryCatch(copula <- setTheta(copula, beta), # try to set parameters
                   error = function(e) NULL)
    if (is.null(tc)) return(-Inf) # in case of failure, return -Inf
    p <- ncol(z) + 1 # number of parameters per margin
    beta.1 <- par[1 + 1:p] # parameters of the first marginal model
    beta.2 <- par[p+1 + 1:p] # parameters of the second marginal model
    beta.3 <- par[2*p+1 + 1:p] # parameters of the third marginal model
    ## Marginal log-likelihood evaluation and computing the
    ## corresponding parametric pseudo-observations
    nmLL.1 <- nmLL(beta.1, x[,1], z, pobs = TRUE)
    nmLL.2 <- nmLL(beta.2, x[,2], z, pobs = TRUE)
    nmLL.3 <- nmLL(beta.3, x[,3], z, pobs = TRUE)
    ## In case of invalid evaluation of the marg. likelihoods, return -Inf
    if (any(is.na(c(nmLL.1$nLL, nmLL.2$nLL, nmLL.3$nLL)))) return(-Inf)
    ## Parametric pseudo-observations as a matrix
    U <- cbind(nmLL.1$U, nmLL.2$U, nmLL.3$U)
    ## -log-likelihood
    -sum(dCopula(u = U, copula = copula, log = TRUE)) +
         nmLL.1$nLL + nmLL.2$nLL + nmLL.3$nLL
}

res.g <- optim(c(ts.g@estimate, ts.math, ts.sci, ts.read), nfLL,
               x = nels, z = z, copula = gumbelCopula(dim = 3),
               method = "BFGS", control = list(maxit = 1e6), hessian = TRUE)
stopifnot(res.g$convergence == 0) # the optimization has converged
-res.g$value # the maximized likelihood

## Maximization of the full conditional likelihood function using
## as initial values the estimates obtained with the two-stage approach
res.f <- optim(c(ts.f@estimate, ts.math, ts.sci, ts.read), nfLL,
               x = nels, z = z, copula = frankCopula(dim = 3),
               method = "BFGS", control = list(maxit = 1e6), hessian = TRUE)
stopifnot(res.f$convergence == 0) # the optimization has converged
-res.f$value # the maximized likelihood

res.n <- optim(c(ts.n.ex@estimate, ts.math, ts.sci, ts.read), nfLL,
               x = nels, z = z, copula = normalCopula(dim = 3),
               method = "BFGS", control = list(maxit = 1e6), hessian = TRUE)
stopifnot(res.n$convergence == 0) # the optimization has converged
-res.n$value # the maximized likelihood

full.n.ex <- res.n$par[1] # parameter estimate of the normal copula
full.math <- res.n$par[1 + 1:p] # param. est. of the 1st marginal model
full.sci  <- res.n$par[p+1 + 1:p] # param. est. of the 2nd marginal model
full.read <- res.n$par[2*p+1 + 1:p] # param. est. of the 3rd marg. model

## Comparison of the estimated copula param. (two-stage vs full likelihood)
c(ts.n.ex@estimate, full.n.ex)

## Comparison of the estimated parameters of the three marginal models
round(cbind(ts.math, full.math, ts.sci, full.sci, ts.read, full.read), 3)

cov.fLL <- solve(res.n$hessian)

sqrt(cov.fLL[1,1])

sqrt(vcov(ts.n.ex))

## Standard errors of marginal parameters (two-stage vs full likelihood)
## Note that, when called on the fitted objects returned by 'glm()',
## 'diag(vcov())' does not provide the variance of the dispersion parameter
full.SE <- sqrt(diag(cov.fLL))
all.SE <- cbind(ts.math   = c(sqrt(diag(vcov(math.glm))), NA),
                full.math = full.SE[1 + 1:p],
                ts.sci    = c(sqrt(diag(vcov(sci.glm))), NA),
                full.sci  = full.SE[p+1 + 1:p],
                ts.read   = c(sqrt(diag(vcov(read.glm))), NA),
                full.read = full.SE[2*p+1 + 1:p])
rownames(all.SE)[p] <- "disp"
round(all.SE, 4)

q.math <- quantile(nels$Math,    probs = 0.1)
q.sci  <- quantile(nels$Science, probs = 0.1)
q.read <- quantile(nels$Reading, probs = 0.1)

probs <- c(0.25, 0.75) # quantile orders of SES
stdnts <- data.frame(Minority = c(0, 1, 0, 1),
                     SES = rep(quantile(nels$SES, probs = probs), each = 2),
                     Female = 0, Public = 1, Size = 0, Urban = 1, Rural = 0)

newz <- model.matrix(~ Minority + SES + Female + Public + Size +
                         Urban + Rural, data = stdnts) # design matrix
## The four marginal probabilities for the math, science and reading score
prob.math <- nmLL(full.math, q.math, newz, pobs = TRUE)$U
prob.sci  <- nmLL(full.sci,  q.sci,  newz, pobs = TRUE)$U
prob.read <- nmLL(full.read, q.read, newz, pobs = TRUE)$U

u <- cbind(prob.math, prob.sci, prob.read)
joint.prob.full <- pCopula(u, copula = normalCopula(full.n.ex, dim = 3))
data.frame(Minority = stdnts[,1] == 1, SES.q.order = rep(probs, each = 2),
           joint.prob = round(joint.prob.full, 6),
           joint.prob.ind = round(pCopula(u, copula = indepCopula(3)), 6))

