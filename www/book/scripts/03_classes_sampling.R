## By Marius Hofert, Ivan Kojadinovic, Martin Maechler, Jun Yan

## R script for Chapter 3 of Elements of Copula Modeling with R

source("00_preliminaries.R")


### 3.1 Elliptical copulas #####################################################

### Construction of elliptical distributions and sampling

## Stepwise generation of a sample from a bivariate t distribution
## via the (general) stochastic representation of elliptical distributions
n <- 1000 # sample size
d <- 2 # dimension
mu <- c(1, 3) # location vector
Sigma <- matrix(c(16, 4,
                   4, 2), ncol = 2, byrow = TRUE) # scale matrix
nu <- 3.5 # degrees of freedom
set.seed(271) # set a seed (for reproducibility)
R <- sqrt(d * rf(n, df1 = d, df2 = nu)) # sample R for a t_nu
A <- t(chol(Sigma)) # Cholesky factor
Z <- matrix(rnorm(n * d), ncol = d) # N_d(0,I_d)
S <- Z/sqrt(rowSums(Z^2)) # uniform distribution on unit sphere (= Z/||Z||)
X <- rep(mu, each = n) + R * t(A %*% t(S)) # X = mu + R A S
plot(S, xlab = quote(S[1]), ylab = quote(S[2]))
plot(t(A %*% t(S)), xlab = quote((AS)[1]), ylab = quote((AS)[2]))
xlim <- range(X[,1], X[,1] - mu[1])
ylim <- range(X[,2], X[,2] - mu[2])
plot(R * t(A %*% t(S)), xlab = quote(R(AS)[1]), ylab = quote(R(AS)[2]),
     xlim = xlim, ylim = ylim)
plot(X, xlab = quote(X[1]), ylab = quote(X[2]),
     xlim = xlim, ylim = ylim)

## 2-point distribution for R and corresponding X
R.2pt <- 1 + rbinom(n, size = 1, prob = 2/3) # prob. 1/3 to be 1, 2/3 to be 2
X.2pt <- rep(mu, each = n) + R.2pt * t(A %*% t(S)) # compute X
## Bounded distribution for R and corresponding X
pR <- function(q) pf(q^2/d, df1 = d, df2 = nu) # df of R for a t_nu
qR <- function(p) sqrt(d * qf(p, df1 = d, df2 = nu)) # quantile function
a <- 1; b <- 3 # [a, b]
R.bdd <- qR(pR(a) + runif(n) * (pR(b) - pR(a))) # sample R on [a,b]
X.bdd <- rep(mu, each = n) + R.bdd * t(A %*% t(S)) # compute X
plot(X.2pt, xlab = quote(X[1]), ylab = quote(X[2]))
plot(X.bdd, xlab = quote(X[1]), ylab = quote(X[2]))


### Normal copula

nc <- normalCopula(iTau(normalCopula(), tau = 0.5))
set.seed(271)
U <- rCopula(1000, copula = nc) # sample from the normal copula
wireframe2(nc, FUN = dCopula, delta = 0.025) # density
contourplot2(nc, FUN = pCopula) # copula
contourplot2(nc, FUN = dCopula, n.grid = 42, cuts = 33, lwd = 1/2) # density
plot(U, xlab = quote(U[1]), ylab = quote(U[2])) # scatter plot


### t copula

nu <- 4 # needs to be an integer here (only) because of pCopula()
tc <- tCopula(iTau(tCopula(, df = nu), tau = 0.5), df = nu)
set.seed(271)
U <- rCopula(1000, copula = tc) # sample from the t copula
wireframe2(tc, FUN = dCopula, delta = 0.025) # density
contourplot2(tc, FUN = pCopula) # copula
contourplot2(tc, FUN = dCopula, n.grid = 42, cuts = 27) # density
plot(U, xlab = quote(U[1]), ylab = quote(U[2])) # scatter plot

## Setup
n <- 1000
d <- 5
nu <- 3.5
rho <- iTau(tCopula(, df = nu), tau = 0.5)
P <- matrix(rho, nrow = d, ncol = d)
diag(P) <- 1
## Method 1: Directly using rmvt() from 'mvtnorm'
library(mvtnorm)
set.seed(271)
X <- rmvt(n, sigma = P, df = nu)
U <- pt(X, df = nu)
## Method 2: Reproducing rmvt()
set.seed(271)
eig <- eigen(P, symmetric = TRUE) # eigenvalue (instead of Cholesky) decomp.
A <- t(eig$vectors %*% (t(eig$vectors) * sqrt(eig$values)))
X.norm <- matrix(rnorm(n * d), nrow = n, byrow = TRUE) %*% A
X.t <- X.norm / sqrt(rchisq(n, nu) / nu)
U.rmvt <- pt(X.t, df = nu)
stopifnot(all.equal(U.rmvt, U))
## Method 3: rCopula()
set.seed(271)
U.cop <- rCopula(n, copula = tCopula(rho, dim = d, df = nu))
stopifnot(all.equal(U.cop, U.rmvt))

X.meta <- cbind(qt(U[,1], df = 2), qt(U[,2], df = 10))
plot(X, xlab = quote(X[1]), ylab = quote(X[2]))
plot(X.meta, xlab = quote(X[1]), ylab = quote(X[2]))


### Marginal inconsistency and its implications

##' @title Random number generator for the exponential power family
##' @param n sample size
##' @param d integer parameter d >= 1
##' @param s numeric parameter s > 0
##' @param gamma numeric parameter gamma > 0
##' @param mu location
##' @param A Cholesky factor of the scale matrix
##' @return Sample from an exponential power family
rExpPow <- function(n, d, s, gamma, mu = rep(0, d), A = diag(d))
{
    R <- (rgamma(n, d / 2 / s) / gamma)^(1 / 2 / s)
    Z <- matrix(rnorm(n * d), ncol = d)
    S <- Z / sqrt(rowSums(Z^2))
    rep(mu, each = n) + R * t(A %*% t(S))
}
## Setup
set.seed(271) # set a seed (for reproducibility)
X.2d <- rExpPow(100000, d = 2, s = 1/2, gamma = 1) # sample for d = 2
X.8d <- rExpPow(100000, d = 8, s = 1/2, gamma = 1) # sample for d = 8
## As the following plot shows, the first univariate margins differ
plot(density(X.2d[,1]), #  est. of 1st univ. marg. of 2-dim. df
     xlim = c(-10, 10), ylim = c(0, 0.35), main = "", xlab = "")
lines(density(X.8d[,1]), col = 2) # est. of 1st univ. marg. of 8-dim. df
legend("topright", bty = "n",
       lty = c(1,1), lwd = c(2, 2), col = 1:2,
       legend = expression("Kernel density estimate of"~X[1]~"(d = 2)",
                           "Kernel density estimate of"~X[1]~"(d = 8)"))
## Empirical marg. prob. trans. followed by quant. trans. to N(0,1)
Z.2d <- qnorm(pobs(X.2d))
Z.8d <- qnorm(pobs(X.8d))
## Kernel est. (on norm. scale) of the copulas of the first bivariate marg.
library(MASS)
dens2d <- kde2d(Z.2d[,1], Z.2d[,2], n = 300)
dens8d <- kde2d(Z.8d[,1], Z.8d[,2], n = 300)
image(dens2d, xlim = c(-2, 2), ylim = c(-2, 2), xlab = quote(Z[1]),
      ylab = quote(Z[2]), col = gray(seq(1, 0, length.out = 100)))
image(dens8d, xlim = c(-2, 2), ylim = c(-2, 2), xlab = quote(Z[1]),
      ylab = quote(Z[2]), col = gray(seq(1, 0, length.out = 100)))


### Grouped normal variance mixture copulas

## Sample from a grouped t copula
d. <- 1:4 # sector dimensions
d <- sum(d.) # total dimension
nu <- rep(4^seq(2, -1), times = d.) # d.o.f. for each sector copula
n <- 1000 # sample size
set.seed(271) # set seed (for reproducibility)
Z <- matrix(rnorm(n * d), ncol = n) # (d, n)-matrix
P <- matrix(0.5, nrow = d, ncol = d) # correlation matrix
diag(P) <- 1 # fix diagonal
A <- t(chol(P)) # Cholesky factor A (s.th. AA^T = P)
Y <- t(A %*% Z) # (n, d) matrix containing n d-vectors following N_d(0, P)
U. <- runif(n) # n-vector of U(0,1) random variates
W <- sapply(nu, function(nu.) 1/qgamma(U., shape = nu./2, rate = nu./2))
X <- sqrt(W) * Y # (n, d)-matrix
U <- sapply(1:d, function(j) pt(X[,j], df = nu[j])) # (n, d)-matrix sample
## Build matrix of colors
cols <- matrix(1, nrow = d, ncol = d) # colors
start <- c(1, cumsum(head(d., n = -1)) + 1) # block start indices
end <- cumsum(d.) # block end indices
for(j in seq_along(d.)) cols[start[j]:end[j], start[j]:end[j]] <- j # colors
diag(cols) <- NA # remove colors corresponding to diagonal entries
splom2(U, pch = ".", pscales = 0, col.mat = cols) # plot sample


### 3.2 Archimedean copulas ####################################################

### Graphs of selected Archimedean generators

copClayton@psi # generator of the Clayton family

t <- seq(0, 2, length.out = 257) # evaluation points
tau <- 0.5 # Kendall's tau
psi. <- cbind(Pi = exp(-t), # Pi generator
              C  = copClayton@psi(t, theta = iTau(claytonCopula(), tau)),
              F  = copFrank@psi  (t, theta = iTau(frankCopula(),   tau)),
              GH = copGumbel@psi (t, theta = iTau(gumbelCopula(),  tau)),
              J  = copJoe@psi    (t, theta = iTau(joeCopula(),     tau)))
plot(t, psi.[,1], type = "l", lwd = 2,
     xlim = range(t), ylim = range(psi.), col = 1, ylab = "",
     xlab = quote(psi(t)~"as a function of t"))
for(j in 2:ncol(psi.)) lines(t, psi.[,j], col = j, lwd = 2)
legend("topright", bty = "n", lty = 1, lwd = 2, col = 1:ncol(psi.),
       legend = c("Independence", "Clayton", "Frank",
                  "Gumbel-Hougaard", "Joe"))


### Sampling from a Gumbel-Hougaard copula

## Setup
n <- 1000 # sample size
d <- 5 # dimension
family <- "Gumbel" # copula family
th <- iTau(archmCopula(family), 0.5) # copula parameter
## Version 1: manually
set.seed(271) # set seed (for reproducibility)
cop <- getAcop(family) # define the Archimedean copula to sample from
V <- cop@V0(n, theta = th) # generate frailties V from F
E <- matrix(-log(runif(n*d)), ncol = d) # sample independent Exp(1)
U.man <- cop@psi(E/V, theta = th) # construct U
## Version 2: via rCopula()
cop <- archmCopula(family, param = th, dim = d) # define the copula to sample
set.seed(271) # set seed (for reproducibility)
U <- rCopula(n, cop) # sample
## Check
stopifnot(all.equal(U.man, U))


### Graphs of |psi^(d)| for the Gumbel-Hougaard generator

t <- 10^seq(-2, 2, by = 0.05) # evaluation points
th <- iTau(gumbelCopula(), 0.5) # corresponding GH parameter theta
d <- c(2, 5, 10, 20, 50) # order of the derivatives
dPsi. <- sapply(d, function(d.) # (length(t), length(d))-mat. of derivatives
    copGumbel@absdPsi(t, theta = th, degree = d.))
plot(t, dPsi.[,1], type = "l", log = "y", lwd = 2,
     xlim = range(t), ylim = range(dPsi.), col = 1, ylab = "",
     xlab = quote(bgroup("|",{psi^(d)}(t),"|")~"as a function of t"))
for(j in 2:length(d)) lines(t, dPsi.[,j], col = j, lwd = 2)
legend("topright", bty = "n", lty = 1, lwd = 2, col = 1:length(d),
       legend = as.expression(lapply(1:length(d), function(j)
           substitute(d==d., list(d.=d[j])))))


### Negative log-likelihood and profile likelihood of a Gumbel-Hougaard copula

## Generate a sample from a Gumbel-Hougaard copula
n <- 100 # sample size
d <- 100 # dimension
tau <- 0.5 # Kendall's tau
family <- "Gumbel" # copula family
th0 <- iTau(archmCopula(family), tau = tau) # true copula parameter
cop <- archmCopula(family, param = th0, dim = d) # define copula
set.seed(271) # set seed (for reproducibility)
U <- rCopula(n, cop) # sample
## Maximum likelihood estimation (MLE)
nLL <- function(th) -loglikCopula(th, u = U, copula = cop) # -log-likelihood
ii <- initOpt(family) # initial interval
iv <- initOpt(family, interval = FALSE, u = U) # initial value
library(bbmle)
mle <- mle2(minuslogl = nLL, optimizer = "optim", method = "L-BFGS-B",
            start = list(th = iv), lower = ii[1], upper = ii[2]) # optimize
stopifnot(mle@details$convergence == 0) # check convergence
th.hat <- mle@coef # estimate
## Confidence intervals (CIs)
CI <- confint(mle, quietly = TRUE)
stopifnot(CI[1] <= th0, th0 <= CI[2]) # check if theta_0 is in the CI
prof <- profile(mle) # profile likelihood (see bbmle's vignette 'mle2')
## nLL plot
th.bds <- iTau(archmCopula(family), tau = c(tau-0.015, tau+0.015)) # bounds
th. <- seq(th.bds[1], th.bds[2], length.out = 101) # x values
nLL. <- unlist(lapply(th., nLL)) # y values
plot(th., nLL., type = "l", xlab = quote(theta), ylab = "")
abline(v = th0, lty = 2)
abline(v = CI[1], lty = 3)
abline(v = CI[2], lty = 3)
abline(v = th.hat, lty = 4)
legend("bottomright", bty = "n", y.intersp = 1.1, lty = 1:4,
       legend = c(quote(- "log-likelihood"),
                  substitute("True"~theta[0]==T0, list(T0 = th0)),
                  "95% CIs", expression("Estimate"~hat(theta)[n])))
mtext(substitute("GH copula,"~~n==N*","~~d==D, list(N = n, D = d)),
      side = 4, line = 0.6, adj = 0, las = 0) # label
## Profile-likelihood plot for theta
plot(prof, xaxt = "n", main = "",
     col.prof = "black", col.minval = "black", col.conf = "black",
     xlabs = expression(theta), # a *vector*
     ylab = quote("Profile likelihood for"~theta * (abs(z)==sqrt(deviance))))
axis(side = 1, at = c(1.96, 1.98, th.hat, th0, 2.02, 2.04),
     labels = c("1.96","1.98", quote(hat(theta)[n]),
                expression(theta[0]), "2.02", "2.04"),
     padj = c(0.13, 0.13, 0.24, 0.35, 0.13, 0.13))


### Outer power Archimedean copulas

## Setup
th.C <- copClayton@iTau(0.3) # Clayton parameter s.t. tau = 0.3
op.C <- opower(copClayton, thetabase = th.C) # define an opC copula family
## Define two opC copulas (tau = 0.5 and 0.8)
th <- sapply(c(0.5, 0.7), op.C@iTau) # choose parameter according to taus
opC  <- onacopulaL(op.C, list(th[1], 1:2)) # define the opC copula
opC. <- onacopulaL(op.C, list(th[2], 1:2))
## Sample
set.seed(271) # set seed (for reproducibility)
U  <- rCopula(1000, copula = opC)
U. <- rCopula(1000, copula = opC.)
plot(U,  xlab = quote(U[1]), ylab = quote(U[2])) # opC sample
plot(U., xlab = quote(U[1]), ylab = quote(U[2])) # opC sample


### Archimedean copulas with d-monotone generators, Liouville copulas

## Setup
n <- 1000 # sample size
d <- 2 # dimension
th <- 1 # Pareto parameter
set.seed(271) # set seed (for reproducibility)
R <- runif(n)^(-1/th) # sample radial part (here: Pareto on [1,Inf))
## Sample from a so-called Pareto-simplex copula
E <- matrix(rexp(n * d), nrow = n, ncol = d) # unit exponentials
S <- E / matrix(rowSums(E), nrow = n, ncol = d) # S uniformly on unit simplex
incBeta <- function(x, a, b) pbeta(x, a, b) * beta(a, b) # incomplete beta
psi <- function(t, th) t^(-1/th) * incBeta(pmin(1,t), a = 1/th, b = d) / th
U <- psi(R * S, th = th)
plot(U, xlab = quote(U[1]), ylab = quote(U[2])) # Pareto-simplex sample

## (Approximately) sample from a so-called Pareto-Liouville copula
alpha <- c(0.5, 2) # Dirichlet parameters
G <- vapply(alpha, function(a) rgamma(n, a), numeric(n)) # Gamma variables
D <- G / rowSums(G) # Dirichlet distribution on unit simplex
U. <- 1 - pobs(R * D) # empirical marginal survival functions
plot(U., xlab = quote(U[1]), ylab = quote(U[2])) # Pareto-Liouville sample


### Nested Archimedean copulas

## Define a nested Gumbel-Hougaard copula
family <- "Gumbel" # copula family
tau <- c(0.2, 0.4, 0.6, 0.8) # Kendall's tau
th <- iTau(archmCopula(family), tau = tau) # corresponding parameters
nlist <- list(th[1], 1, list(list(th[2], 2:3), # NAC structure
                             list(th[3], 4, list(list(th[4], 5:7)))))
NAC <- onacopulaL(family, nacList = nlist) # NAC copula
## Sample
set.seed(271) # set seed (for reproducibility)
U <- rCopula(1000, copula = NAC) # sample
## Build matrix of colors
cols <- matrix(1, nrow = 7, ncol = 7)
cols[2, 3] <- 2
cols[4, 5:7] <- 3
cols[5:7, 5:7] <- 4
cols[lower.tri(cols)] <- t(cols)[lower.tri(cols)]
diag(cols) <- NA
splom2(U, pch = ".", pscales = 0, col.mat = cols) # plot sample


### 3.3 Extreme-value copulas ##################################################

### Exchangeable extreme-value copulas

## Parameter values corresponding to a Kendall's tau of 0.25, 0.5 and 0.75
th.25 <- iTau(gumbelCopula(), tau = 0.25)
th.50 <- iTau(gumbelCopula(), tau = 0.50)
th.75 <- iTau(gumbelCopula(), tau = 0.75)
## Graphs of the corresponding Pickands dependence functions A
curve(A(gumbelCopula(th.25), w = x), from = 0, to = 1, ylim = c(0.5 ,1),
      xlab = "t", ylab = expression({A[theta]^{GH}}(t)), col = 1, lwd = 2)
curve(A(gumbelCopula(th.50), w = x), add = TRUE, col = 2, lwd = 2)
curve(A(gumbelCopula(th.75), w = x), add = TRUE, col = 3, lwd = 2)
## Every bivariate Pickands dependence function is convex
## and its graph is necessarily in the "triangle" plotted below
lines(c(0, 1), c(1, 1), lty = 2)
lines(c(0, 0.5, 1), c(1, 0.5, 1), lty = 2)
legend("bottomright", bty = "n", lwd = 2, col = 1:3,
       legend = expression(tau == 0.25, tau == 0.5, tau == 0.75),
       inset = 0.02, y.intersp = 1.2)

## Draw the Pickands dependence function of the independence copula
## which is the same as "lines(c(0, 1), c(1, 1), lty = 2)"
curve(A(indepCopula(), w = x), from = 0, to = 1, ylim = c(0.5, 1),
      xlab = "t", ylab = quote({A[theta]}(t)), lty = 2)
## Add the constraint related to the monotonicity copula
lines(c(0, 0.5, 1), c(1, 0.5, 1), lty = 2)
## Set Kendall's tau and color
tau <- 0.25
col <- adjustcolor(1, alpha.f = 0.25)
## And the utility A() curve plot function
curveA <- function(x, Cop) # using  global (tau, col)
  curve(A(Cop(iTau(Cop(), tau)), w = x), add = TRUE, col = col, lwd = 2)

curveA(x, galambosCopula)
curveA(x, gumbelCopula)
curveA(x, huslerReissCopula)
curveA(x, tevCopula)
if(tau < 0.4184) curveA(x, tawnCopula) # no tawnCopula for larger tau!


### 3.4 Selected copula transformations and constructions ######################

### 3.4.1 Rotated copulas ######################################################

### A rotated bivariate Clayton copula

## The vector r is represented by a vector of logicals
rc <- rotCopula(claytonCopula(4), flip = c(TRUE, FALSE))
wireframe2(rc, FUN = pCopula, # wireframe plot (copula)
           draw.4.pCoplines = FALSE)
wireframe2(rc, FUN = dCopula, delta = 0.025) # wireframe plot (density)
contourplot2(rc, FUN = pCopula, n.grid = 64) # contour plot (copula)
contourplot2(rc, FUN = dCopula, n.grid = 64, cuts = 30,
             pretty = FALSE, lwd = 1/2) # contour plot (density)


### A rotated four-dimensional Frank copula

## The logical representing vector r
flip <- c(TRUE, FALSE, TRUE, FALSE)
rf <- rotCopula(frankCopula(10, dim = 4), flip = flip)
n <- 1000
set.seed(2016)
U <- rCopula(n, copula = rf)
set.seed(2016)
V <- rCopula(n, frankCopula(10, dim = 4))
## "Flip" the relevant components
V[, flip] <- 1 - V[, flip]
stopifnot(all.equal(U, V)) # check

splom2(U, cex = 0.3, col.mat = "black")


### 3.4.2 Khoudraji's device ###################################################

### Non-exchangeable Khoudraji-Clayton copulas

(kc <- khoudrajiCopula(copula2 = claytonCopula(6), shapes = c(0.2, 0.95)))

n <- 5000
plot(rCopula(n, copula = kc),
     cex = 0.5, xlab = quote(U[1]), ylab = quote(U[2]))
plot(rCopula(n, copula = setTheta(kc, value = c(6, 0.4, 0.95))),
     cex = 0.5, xlab = quote(U[1]), ylab = quote(U[2]))
plot(rCopula(n, copula = setTheta(kc, value = c(6, 0.6, 0.95))),
     cex = 0.5, xlab = quote(U[1]), ylab = quote(U[2]))
plot(rCopula(n, copula = setTheta(kc, value = c(6, 0.8, 0.95))),
     cex = 0.5, xlab = quote(U[1]), ylab = quote(U[2]))

plot(rCopula(n, copula = setTheta(kc, value = c(6, 0.95, 0.6))),
     cex = 0.5, xlab = quote(U[1]), ylab = quote(U[2]))

plot(rCopula(n, copula = khoudrajiCopula(copula1 = claytonCopula(6),
                                         shapes = c(0.6, 0.95))),
     cex = 0.5, xlab = quote(U[1]), ylab = quote(U[2]))

plot(rCopula(n, copula = khoudrajiCopula(copula1 = claytonCopula(6),
                                         shapes = c(1 - 0.6, 1 - 0.95))),
     cex = 0.5, xlab = quote(U[1]), ylab = quote(U[2]))


### Non-exchangeable Khoudraji-Gumbel-Hougaard-Clayton copulas

## The setup
s <- c(0.6, 0.95)
copula1 <- gumbelCopula
copula2 <- claytonCopula
## A utility function to obtain the parameter values of C_1 and C_2
param <- function(tau) c(iTau(copula1(), tau), iTau(copula2(), tau))
## The corresponding Khoudraji-Gumbel-Hougaard-Clayton copula
(kho <- khoudrajiCopula(copula1 = copula1(param(0.65)[1]),
                        copula2 = copula2(param(0.65)[2]),
                        shapes = s))

n <- 5000
U <- rCopula(n, copula = kho)
plot(U, cex = 0.5, xlab = quote(U[1]), ylab = quote(U[2]))
V <- rCopula(n, copula = setTheta(kho, value = c(param(0.8), s)))
plot(V, cex = 0.5, xlab = quote(V[1]), ylab = quote(V[2]))
W <- rCopula(n, copula = setTheta(kho, value = c(param(0.95), s)))
plot(W, cex = 0.5, xlab = quote(W[1]), ylab = quote(W[2]))

c(cor(U, method = "kendall")[1,2], cor(V, method = "kendall")[1,2],
  cor(W, method = "kendall")[1,2])


### A non-exchangeable extreme-value family

kg <- khoudrajiCopula(copula2 = gumbelCopula(4), shapes = c(0.2, 0.95))
curve(A(kg, w = x), from = 0, to = 1, ylim = c(0.5, 1),
      xlab = "t", ylab = expression({A[theta]^{KGH}}(t)), col = 1, lwd = 2)
curve(A(setTheta(kg, value = c(4, 0.4, 0.95)), w = x), add = TRUE, lwd = 2,
      col = 2)
curve(A(setTheta(kg, value = c(4, 0.6, 0.95)), w = x), add = TRUE, lwd = 2,
      col = 3)
curve(A(setTheta(kg, value = c(4, 0.8, 0.95)), w = x), add = TRUE, lwd = 2,
      col = 4)
lines(c(0, 1),      c(1, 1),      lty = 2)
lines(c(0, 0.5, 1), c(1, 0.5, 1), lty = 2)
legend("bottomright", bty = "n", lwd = 2,, col = 1:4,
       legend = expression(s[1] == 0.2, s[1] == 0.4,
                           s[1] == 0.6, s[1] == 0.8))


### Higher-dimensional Khoudraji copulas

kgsc <- khoudrajiCopula(copula1 = gumbelCopula(2, dim=3),
                        copula2 = rotCopula(claytonCopula(6, dim=3)),
                        shapes = c(.6, 0.7, 0.95))
## Random points in the unit hypercube where to evaluate the density
set.seed(42)
v <- matrix(runif(15), 5, 3)
dCopula(v, copula = kgsc)

kgn <- khoudrajiCopula(copula1 = gumbelCopula(2, dim=3),
                       copula2 = normalCopula(0.9, dim=3),
                       shapes = c(.6, 0.7, 0.95))
try( dCopula(v, copula = kgn) ) # not implemented


### 3.4.3 Mixtures of copulas ##################################################

### A mixture of Clayton and Gumbel-Hougaard copulas

cc <- claytonCopula(iTau(claytonCopula(), tau = 0.75)) # the first component
gc <- gumbelCopula(iTau(gumbelCopula(),   tau = 0.75)) # the second component
weights <- c(1/3, 2/3) # the corresponding weights
(mcg <- mixCopula(list(cc, gc), w = weights)) # the mixture copula

stopifnot(all.equal(rho(mcg), weights[1] * rho(cc) + weights[2] * rho(gc)))
lambda(mcg)

stopifnot(all.equal(lambda(mcg),
                    weights[1] * lambda(cc) + weights[2] * lambda(gc)))

set.seed(127)
U <- rCopula(1000, copula = mcg) # sample from the mixture
wireframe2(mcg, FUN = dCopula, delta = 0.025) # density
contourplot2(mcg, FUN = pCopula) # copula
contourplot2(mcg, FUN = dCopula, cuts = 32, # density
             n.grid = 50, pretty = FALSE,
             col = adjustcolor(1, 1/3), alpha.regions = 3/4)
plot(U, xlab = quote(U[1]), ylab = quote(U[2])) # scatter plot

