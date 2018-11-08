## By Marius Hofert, Ivan Kojadinovic, Martin Maechler, Jun Yan

## R script for Chapter 4 of Elements of Copula Modeling with R


library(copula)


### 4.1 Estimation under a parametric assumption on the copula #################

### 4.1.1 Parametrically estimated margins #####################################

### Estimation of copula parameters via the MLE

## The "unknown" copula (a 2-dim. Clayton copula with parameter 3)
cc <- claytonCopula(3)
## The "unknown" distribution (N(0,1), Exp(1) margins)
mcc <- mvdc(cc, margins = c("norm", "exp"),
            paramMargins = list(list(mean = 0, sd = 1),
                                list(rate = 1)))
## Generate the "observed" sample
set.seed(712)
n <- 1000
X <- rMvdc(n, mvdc = mcc)
## The function fitMvdc() estimates all the parameters of the mvdc object
## mcc (whose parameter values are not used). Starting values need to be
## provided.
start <- c(mu0 = mean(X[,1]), sig0 = sd(X[,1]), lam0 = 1 / mean(X[,2]),
           th0 = 2)
(mle <- fitMvdc(X, mvdc = mcc, start = start))

summary(mle)


### Estimation of copula parameters via the IFME

## Parametric pseudo-observations obtained from X by marginal MLE
U <- cbind(pnorm(X[,1], mean = mean(X[,1]),
                 sd = sqrt((n - 1) / n) * sd(X[,1])),
           pexp(X[,2], rate = 1 / mean(X[,2])))
ifme <- fitCopula(claytonCopula(), data = U, method = "ml")
summary(ifme)

optimMeth(claytonCopula(), method = "ml", dim = 2)


### 4.1.2 Nonparametrically estimated margins #################################

### Pseudo-observations of daily log-returns

data(rdj) # 'head(rdj)' for looking at the first six observations
splom2(rdj[,2:4], cex = 0.4, col.mat = adjustcolor("black", 0.5))

U <- pobs(rdj[,2:4])
splom2(U, cex = 0.4, col.mat = "black")


### Estimation of copula parameters via the method of moments based on Kendall's tau

## The "unknown" copula (a 2-dim. Gumbel-Hougaard copula with parameter 3)
gc <- gumbelCopula(3)
## The "unknown" distribution (N(0,1) margins)
mgc <- mvdc(gc, margins = c("norm", "norm"),
            paramMargins = list(list(mean = 0, sd = 1),
                                list(mean = 0, sd = 1)))
## Generate the "observed" sample
set.seed(49)
X <- rMvdc(1000, mvdc = mgc)
## The sample version of Kendall's tau
tau.n <- cor(X[,1], X[,2], method = "kendall")
## The corresponding copula parameter estimate
(itau <- iTau(gc, tau = tau.n))

stopifnot(all.equal(itau, 1 / (1 - tau.n))) # the same
## The same but with a standard error
summary(fitCopula(gumbelCopula(), data = pobs(X), method = "itau"))


### Estimation of copula parameters via the method of moments based on Spearman's rho

## The "unknown" copula (a 2-dim. normal copula with parameter 0.5)
nc <- normalCopula(0.5)
## Generate the "observed" sample
set.seed(314)
X <- rCopula(1000, nc)
## The sample estimate of Spearman's rho
rho.n <- cor(X[,1], X[,2], method = "spearman")
## The corresponding copula parameter estimate
(irho <- iRho(nc, rho = rho.n))

stopifnot(all.equal(irho, 2 * sin(pi * rho.n / 6))) # the same
## The same but with a standard error
summary(fitCopula(normalCopula(), data = pobs(X), method = "irho"))


### Application of the estimation via the method of moments based on Kendall's tau

data(danube, package = "lcopula") # already pseudo-observations
U <- as.matrix(danube)
plot(U, xlab = "Donau", ylab = "Inn")

summary(fitCopula(gumbelCopula(), data = U, method = "itau"))

summary(fitCopula(plackettCopula(), data = U, method = "itau"))

summary(fitCopula(normalCopula(), data = U, method = "itau"))


### Estimation of copula parameters via the MPLE

## The "unknown" copula (a 2-dim. Frank copula with parameter 3)
fc <- frankCopula(3)
## Generate the "observed" sample
set.seed(271)
U <- rCopula(1000, fc)
## Compute the MPLE and its standard error
summary(fitCopula(frankCopula(), data = pobs(U), method = "mpl"))

optimMeth(frankCopula(), method = "mpl", dim = 2)


### Application of the estimation via the MPLE

U <- pobs(rdj[,2:4]) # compute the pseudo-observations
## MPLE for the normal copula
summary(fitCopula(normalCopula(dim = 3, dispstr = "un"), data = U))

## MPLE for the t copula
summary(fitCopula(tCopula(dim = 3, dispstr = "un"), data = U))


### 4.1.3 Estimators of elliptical copula parameters ###########################

### Estimation of normal copula parameters via method-of-moments

data(SMI.12) # load the SMI constituent data
library(qrmtools)
X <- returns(SMI.12) # compute log-returns
U <- pobs(X) # compute pseudo-observations
d <- ncol(U) # 20 dimensions

f.irho <- fitCopula(normalCopula(dim = d, dispstr = "un"), data = U,
                    method = "irho")
f.itau <- fitCopula(normalCopula(dim = d, dispstr = "un"), data = U,
                    method = "itau")

P.irho <- p2P(coef(f.irho), d = d)
P.itau <- p2P(coef(f.itau), d = d)


### Estimation of t copula parameters using the method of Mashal and Zeevi

fit <- fitCopula(tCopula(dim = d, dispstr = "un"), data = U,
                 method = "itau.mpl")


### Linear correlation vs Kendall's tau for estimating $t$ distributions

##' @title Compute (repeated) parameter estimates for the bivariate t
##' @param nu degrees of freedom parameter
##' @param rho correlation parameter
##' @param B number of replications
##' @param n sample size
##' @return (B, 2)-matrix containing the B estimates computed
##'         via Pearson's rho and Kendall's tau
estimates <- function(nu, rho, B, n)
{
    ## Generate data (B-list of (n, 2)-matrices)
    tc <- tCopula(rho, df = nu) # define the corresponding t copula
    X <- lapply(1:B, function(i) qt(rCopula(n, copula = tc), df = nu))
    ## For each of the B data sets, estimate the correlation parameter of the
    ## t distribution via Pearson's rho and Kendall's tau.
    pea <- vapply(X, function(x) cor(x[,1], x[,2]), numeric(1))
    ken <- iTau(tCopula(df = nu),
                tau = sapply(X, function(x) cor(x, method = "kendall")[2,1]))
    ## Return
    cbind(Pearson = pea, Kendall = ken) # (B, 2)-matrix
}

set.seed(271) # for reproducibility
nu <- 3 # degrees of freedom
rho <- 0.5 # correlation parameter
r <- estimates(nu, rho = rho, B = 3000, n = 90) # (B, 2)-matrix
varP <- var(r[,"Pearson"]) # variance of sample linear correlation
varK <- var(r[,"Kendall"]) # variance of inverting Kendall's tau
VRF <- varP / varK # variance reduction factor
PIM <- (varP - varK) / varP * 100 # % improvement
boxplot(r, names = c("Sample linear correlation", "Inverting Kendall's tau"),
        ylab = substitute("Estimates of"~rho~"of a"
                          ~t[nu.]~"distribution with"~rho==rho.,
                          list(nu. = nu, rho. = rho)))
mtext(substitute("VRF (% improvement):"~~v~"("*i*"%)",
                 list(v = round(VRF, 2), i = round(PIM))),
      side = 4, line = 1, adj = 0, las = 0)

## Estimate variance reduction factor for various nu and rho
nu. <- 2 + 2^seq(-2, 7, by = 0.5) # degrees of freedom to consider
rho. <- c(-0.99, -0.8, 0, 0.8, 0.99) # correlations to consider
res <- lapply(nu., function(nu..) lapply(rho., function(rho..)
              estimates(nu.., rho = rho.., B = 3000, n = 90)))
vars <- lapply(res, function(r) lapply(r, function(r.)
               var(r.[,"Pearson"])/var(r.[,"Kendall"])))
VRF. <- matrix(unlist(vars), nrow = length(rho.), ncol = length(nu.),
               dimnames = list(rho = rho., nu = nu.))
ylim <- range(VRF.)
lrho <- length(rho.)
plot(nu., VRF.[1,], type = "l", log = "xy", ylim = ylim, col = 1,
     xlab = "Degrees of freedom",
     ylab = "VRF (correlation over inverting Kendall's tau)")
for(k in 2:lrho) lines(nu., VRF.[k,], col = k)
abline(h = 1, lty = 2)
legend("topright", bty = "n", lty = rep(1, lrho), col = 1:lrho,
       legend = as.expression(lapply(1:lrho, function(k)
           substitute(rho == rho., list(rho. = rho.[k])))))


### 4.1.5 Estimation of copula models with partly fixed parameters #############

### Estimation of elliptical copulas with partly fixed parameters

## The "unknown" copula (a 3-dim. normal copula)
nc  <- normalCopula(param = c(0.6, 0.3, 0.2), dim = 3, dispstr = "un")
## Generate the "observed" sample and compute corresponding pobs
set.seed(819)
U <- pobs(rCopula(1000, nc))
## A trivariate normal copula whose first parameter is fixed to 0.6
(ncf <- normalCopula(param = fixParam(c(0.6, NA_real_, NA_real_),
                                      c(TRUE, FALSE, FALSE)),
                     dim = 3, dispstr = "un"))

fitCopula(ncf, data = U) # MPLE

fitCopula(ncf, data = U, method = "itau")

fitCopula(ncf, data = U, method = "irho")

fixedParam(nc) <- c(TRUE, FALSE, FALSE)
nc

## The "unknown" copula is a 3-dim. t copula (with 4 d.o.f. by default)
tc  <- tCopula(param = c(0.6, 0.3, 0.2), dim = 3, dispstr = "un")
## Generate the "observed" sample and compute corresponding pobs
set.seed(314)
U <- pobs(rCopula(1000, tc))
## A trivariate t copula whose first two parameters are fixed to 0.6 and 0.3
(tcf <- tCopula(param = fixParam(c(0.6, 0.3, NA_real_),
                                 c(TRUE, TRUE, FALSE)),
                dim = 3, dispstr = "un"))

fitCopula(tcf, data = U) # MPLE

fitCopula(tcf, data = U, method = "itau.mpl")

## A trivariate t copula whose first correlation parameter is fixed to 0.6
## and whose number of degrees of freedom is fixed to 4 (default value)
(tcf2 <- tCopula(param = fixParam(c(0.6, NA_real_, NA_real_),
                                  c(TRUE, FALSE, FALSE)),
                 dim = 3, dispstr = "un", df.fixed = TRUE))

fitCopula(tcf2, data = U) # MPLE

fitCopula(tcf2, data = U, method = "itau")


### Estimation of Khoudraji-Clayton copulas with partly fixed parameters

## The "unknown" copula (a 2-dim. Khoudraji-Clayton copula)
kc <- khoudrajiCopula(copula2 = claytonCopula(6), shapes = c(0.4, 1))
set.seed(1307)
U <- pobs(rCopula(1000, kc))
## The default optimization method for fitting bivariate copulas
## constructed with Khoudraji's device by MPLE
optimMeth(khoudrajiCopula(), method = "mpl", dim = 2)

try(fitCopula(khoudrajiCopula(copula2 = claytonCopula()),
              start = c(1.1, 0.5, 0.5), data = U))

fitCopula(khoudrajiCopula(copula2 = claytonCopula()),
          start = c(1.1, 0.5, 0.5), data = U,
          optim.method = "Nelder-Mead")

kcf <- khoudrajiCopula(copula2 = claytonCopula(),
                       shapes = fixParam(c(NA_real_, 1),
                                         c(FALSE, TRUE)))
fitCopula(kcf, start = c(1.1, 0.5), data = U)


### 4.2 Nonparametric estimation of the copula ################################

### 4.2.1 The empirical copula #################################################

### Nonparametric estimation by the empirical copula

## The "unknown" copula (a 3-dim. Clayton copula with parameter 3)
d <- 3
cc <- claytonCopula(3, dim = d)
## Generate a sample from the copula, which will be transformed
## to pseudo-observations in 'C.n()'
n <- 1000
set.seed(65)
U <- rCopula(n, copula = cc)
## Generate random points where to evaluate the empirical copula
v <- matrix(runif(n * d), nrow = n, ncol = d)
ec <- C.n(v, X = U)
## Compare with the true copula; increase n to decrease the error
true <- pCopula(v, copula = cc)
round(mean(abs(true - ec) / true) * 100, 2) # mean relative error (in %)


### The empirical beta and checkerboard copulas

gc <- gumbelCopula(4, dim = 3) # the 'unknown' copula
##' @title Mean relative error in % (for empirical copula, empirical beta
##'        copula and empirical checkerboard copula)
##' @param n sample size
##' @param cop copula object
##' @return mean relative errors in %
compareEmpCops <- function(n, cop)
{
    d <- dim(cop)
    U <- rCopula(n, copula = cop) # a sample from the true copula
    v <- matrix(runif(n * d), nrow = n, ncol = d) # random evaluation points
    ec    <- C.n(v, X = U) # the empirical copula values
    beta  <- C.n(v, X = U, smoothing = "beta") # the emp. beta cop. values
    check <- C.n(v, X = U, smoothing = "checkerboard") # emp. check. cop. val
    true <- pCopula(v, copula = cop) # the true copula values
    c(ec    = mean(abs(true - ec) / true),
      beta  = mean(abs(true - beta) / true),
      check = mean(abs(true - check) / true)) * 100 # mean rel. error in %
}

set.seed(2013)
round(rowMeans(replicate(100, compareEmpCops(30, cop = gc))), 2)

set.seed(2013)
round(rowMeans(replicate(100, compareEmpCops(300, cop = gc))), 2)

set.seed(2008)
U <- rCopula(30, copula = gc) # a sample from the true copula
m <- 100 # number of evaluation points
v <- runif(m) # random points where to evaluate the margins of the estimators
w1 <- cbind(v, 1, 1) # evaluations points margin 1
w2 <- cbind(1, v, 1) # evaluations points margin 2
w3 <- cbind(1, 1, v) # evaluations points margin 3
## Checks
stopifnot(all.equal(C.n(w1, X = U, smoothing = "beta"), v))
stopifnot(all.equal(C.n(w2, X = U, smoothing = "beta"), v))
stopifnot(all.equal(C.n(w3, X = U, smoothing = "beta"), v))
stopifnot(all.equal(C.n(w1, X = U, smoothing = "checkerboard"), v))
stopifnot(all.equal(C.n(w2, X = U, smoothing = "checkerboard"), v))
stopifnot(all.equal(C.n(w3, X = U, smoothing = "checkerboard"), v))


### 4.2.2 Under extreme-value dependence #######################################

### Nonparametric estimation of the Pickands dependence function

## The "unknown" copula (a 2-dim. extreme-value copula)
kg <- khoudrajiCopula(copula1 = indepCopula(),
                      copula2 = gumbelCopula(3),
                      shapes = c(0.6, 0.95))
## Generate a sample from this copula transformed to pseudo-observations
set.seed(172)
U <- pobs(rCopula(100, copula = kg))

## Graphs of the Pickands dependence function A and of its two estimates
curve(An.biv(U, x, estimator = "Pickands"), from = 0, to = 1, col = 1,
      lwd = 2, ylim = c(0.5, 1), xlab = "t", ylab = "A(t)")
curve(An.biv(U, x, estimator = "CFG"), 0, 1, add = TRUE, lwd = 2, col = 2)
curve(A(kg, w = x), 0, 1, add = TRUE, lwd = 2, col = 3)
lines(c(0, 0.5, 1), c(1, 0.5, 1), lty = 2)
lines(c(0, 1),      c(1, 1),      lty = 2)
legend("bottomright", bty = "n", lwd = 2, col = 1:3,
       legend = expression({A[list(n,c)]^{P}}(t), {A[list(n,c)]^{CFG}}(t),
                           {A[bold(theta)]^{KGH}}(t)),
       inset = 0.02, y.intersp = 1.2)

