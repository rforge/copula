## By Marius Hofert, Ivan Kojadinovic, Martin Maechler, Jun Yan

## R script for Chapter 5 of Elements of Copula Modeling with R


### 5.1 Basic graphical diagnostics ############################################

### Pseudo-observations and normal scores

data(danube, package = "lcopula")
U <- as.matrix(danube)
plot(U, xlab = "Donau", ylab = "Inn")
plot(qnorm(U), xlab = "Donau", ylab = "Inn")


### Comparing (non-)parametric estimates of the copula

## Fit a Gumbel-Hougaard copula and compute the contours
fg <- fitCopula(gumbelCopula(), data = U)
cpG <- contourplot2(fg@copula, FUN = pCopula, region = FALSE,
                    key = list(corner = c(0.04, 0.04),
                               lines = list(col = 1:2, lwd = 2),
                               text = list(c("Fitted Gumbel-Hougaard copula",
                                             "Empirical copula"))))
## Fit a Joe copula and compute the contours
fj <- fitCopula(joeCopula(), data = U)
cpJ <- contourplot2(fj@copula, FUN = pCopula, region = FALSE,
                    key = list(corner = c(0.04, 0.04),
                               lines = list(col = 1:2, lwd = 2),
                               text = list(c("Fitted Joe copula",
                                             "Empirical copula"))))
## Compute the contours of the empirical copula
u <- seq(0, 1, length.out = 16)
grid <- as.matrix(expand.grid(u1 = u, u2 = u))
val <- cbind(grid, z = C.n(grid, X = U))
cpCn <- contourplot2(val, region = FALSE, labels = FALSE, col = 2)
## Plots (lattice objects)
library(latticeExtra)
cpG + cpCn
cpJ + cpCn

(fk <- fitCopula(khoudrajiCopula(copula2 = gumbelCopula()), data = U,
                 start = c(1.1, 0.5, 0.5), optim.method = "Nelder-Mead"))

(fk2 <- fitCopula(khoudrajiCopula(copula2 = gumbelCopula(),
                                  shapes = fixParam(c(NA_real_, 1),
                                                    c(FALSE, TRUE))),
                  data = U, start = c(1.1, 0.5)))

cpK <- contourplot2(fk2@copula, FUN = pCopula, region = FALSE,
                    key = list(corner = c(0.04, 0.04),
                               lines = list(col = 1:2, lty = 1),
                               text = list(c("Fitted Khoudraji-Gumbel-Hougaard copula",
                                             "Empirical copula"))))
cpK + cpCn


### Comparing (non-)parametric estimates of the Pickands dependence function

curve(A(fg@copula, x), from = 0, to = 1, ylim = c(0.5, 1),
      lwd = 2, xlab = "t", ylab = "A(t)", col = 1) # parametric
curve(A(fk2@copula, x), 0,1, lwd = 2, add = TRUE, col = 2) # parametric
curve(An.biv(U, x),     0,1, lwd = 2, add = TRUE, col = 3) # non-parametric
lines(c(0, 0.5, 1), c(1, 0.5, 1), lty = 2)
lines(c(0, 1),      c(1, 1),      lty = 2)
legend("bottomright", bty = "n", lwd = 2, col=1:3,
       legend = expression({A[theta[n]]^{GH}}(t),
                           {A[bold(theta)[n]]^{KGH}}(t),
                           {A[list(n,c)]^{CFG}}(t)),
       inset = 0.02, y.intersp = 1.2)


### 5.2 Hypothesis tests #######################################################

### 5.2.1 Tests of independence ################################################

### Test of uncorrelatedness

data(danube, package = "lcopula")
cor.test(danube[,1], danube[,2], method = "kendall")


### A fallacy of a test of uncorrelatedness

set.seed(1515)
x <- rnorm(200)
y <- x^2
cor.test(x, y, method = "kendall")


### Test of independence based on S_n^Pi

n <- 100
d <- 3
set.seed(1969)
U <- rCopula(n, frankCopula(2, dim = d))

dist <- indepTestSim(n, p = d, verbose = FALSE)
indepTest(U, d = dist)


### 5.2.2 Tests of exchangeability #############################################

### Test of exchangeability based on S_n^{ex_C}

set.seed(1453)
exchTest(as.matrix(danube))


### Test of exchangeability based on S_n^{ex_A}

withTime <- function(expr, ...)
{
    st <- system.time(r <- expr, ...)
    list(value = r,  sys.time = st)
}

set.seed(1492)
withTime(exchEVTest(as.matrix(danube)))


### 5.2.3 A test of radial symmetry ############################################

### Test of radial symmetry based on S_n^sym

set.seed(1453)
withTime(radSymTest(as.matrix(danube)))

data(rdj)
Xrdj <- as.matrix(rdj[,-1]) # omitting component 'Date'
set.seed(1389)
withTime(radSymTest(Xrdj))


### 5.2.4 Tests of extreme-value dependence ####################################

### Test of extreme-value dependence based on S_n^{ev_K}

set.seed(1805)
evTestK(as.matrix(danube))


### Test of extreme-value dependence based on S_n^{ev_A}

set.seed(1815)
withTime(evTestA(as.matrix(danube)))


### Test of extreme-value dependence based on S_n^{ev_C}

set.seed(1905)
withTime(evTestC(Xrdj))


### 5.2.5 Goodness-of-fit tests ################################################

### Parametric bootstrap-based tests

set.seed(1598)
withTime(gofCopula(claytonCopula(dim = 3), x = Xrdj, optim.method="BFGS"))


### Multiplier goodness-of-fit tests

set.seed(1610)
withTime(gofCopula(claytonCopula(dim = 3), x = Xrdj, simulation = "mult"))

set.seed(1685)
gofCopula(normalCopula(dim = 3, dispstr = "un"), x = Xrdj,
          simulation = "mult")

set.seed(1792)
gofCopula(tCopula(dim = 3, dispstr = "un", df.fixed = TRUE, df = 10),
          x = Xrdj, simulation = "mult")


### Empirical levels of the multiplier goodness-of-fit test for the Joe family

theta <- iTau(joeCopula(), tau = 0.5) # Joe copula parameter
##' @title P-value of multiplier goodness-of-fit test on data
##'        generated under the null hypothesis
##' @param n sample size
##' @param theta Joe copula parameter
##' @return p-value of the multiplier goodness-of-fit test
pvalMult <- function(n, theta)
{
    U <- rCopula(n, copula = joeCopula(theta))
    gofCopula(joeCopula(), x = pobs(U), simulation = "mult",
              optim.method = "BFGS")$p.value
}

set.seed(1940)
pv <- withTime(replicate(1000, pvalMult(n = 100, theta = theta)))

pv$sys.time # the run time

alpha <- c(0.01, 0.05, 0.1) # nominal levels
rbind(nom.level = alpha, emp.level = ecdf(pv$value)(alpha)) # levels


### Goodness-of-fit tests based on method-of-moments estimation

Xd <- as.matrix(danube)
set.seed(1613)
gofCopula(gumbelCopula(), x = Xd, estim.method = "itau")

set.seed(1914)
gofCopula(gumbelCopula(), x = Xd, estim.method = "itau",
          simulation = "mult")

set.seed(1848)
gofCopula(joeCopula(), x = Xd, estim.method = "itau")


### Goodness of fit of the Khoudraji-Gumbel-Hougaard family

set.seed(1969) ## Parametric bootstrap-based test
withTime(
    gofCopula(khoudrajiCopula(copula2 = gumbelCopula(),
                          shapes = fixParam(c(NA_real_, 1), c(FALSE, TRUE))),
              start = c(1.1, 0.5), x = Xd, optim.method = "Nelder-Mead")
)

set.seed(1969) ## Multiplier-based test
withTime(
    gofCopula(khoudrajiCopula(copula2 = gumbelCopula(),
                          shapes = fixParam(c(NA_real_, 1), c(FALSE, TRUE))),
              start = c(1.1, 0.5), x = Xd, optim.method = "Nelder-Mead",
              simulation = "mult")
)


### 5.2.6 A mixture of graphical and formal goodness-of-fit tests ##############

### A mixture of graphical and formal goodness-of-fit tests

## Load the data, compute the log-returns and the pseudo-observations
data(SMI.12) # load the SMI constituent data
library(qrmtools)
X <- returns(SMI.12) # compute log-returns
U <- pobs(X) # compute pseudo-observations
d <- ncol(U) # 20 dimensions

fit <- fitCopula(tCopula(dim = d, dispstr = "un"), data = U,
                 method = "itau.mpl")
## Extract parameter estimates
len <- length(fit@estimate)
stopifnot(len == d*(d-1)/2 + 1) # sanity check
p <- fit@estimate[seq_len(len-1)] # correlations
## Note: The estimated correlation matrix can be obtained via p2P(p, d = d)
nu <- fit@estimate[len] # degrees of freedom nu

## Define the H_0 copula
cop.t <- ellipCopula("t", df = nu, param = p, dim = d, dispstr = "un")
## Build the array of pairwise H_0-transformed data columns
cu.u.t <- pairwiseCcop(U, cop.t, df = nu)
## Compute pairwise matrix (d by d) of p-values and corresponding colors
set.seed(1389) # for reproducibility
pw.indep.t <- pairwiseIndepTest(cu.u.t, N = 256, verbose = FALSE)
p.val.t <- pviTest(pw.indep.t) # matrix of p-values
## Pairwise Rosenblatt plot (test for the fitted t copula)
title <- list("Pairwise Rosenblatt transformed pseudo-observations",
              quote(bold("to test")~~italic(H[0]:C)~~"is"~~italic(t)))
cols <- gray(seq(0.1, 1, length.out = 6))
pairsRosenblatt(cu.u.t, pvalueMat = p.val.t, method = "QQchisq", pch=".",
                xaxt = "n", yaxt = "n",
                colList = pairsColList(p.val.t, bucketCols = cols,
                                       BWcutoff = 0),
                main.centered = TRUE, main = title, line.main = c(2, -0.8))


### 5.3 Model selection ########################################################

### Cross-validation for the danube data set

Xdan <- as.matrix(danube)
withTime(xvCopula(joeCopula(), x = Xdan))

withTime(xvCopula(gumbelCopula(), x = Xdan))

withTime(
    xvCopula(khoudrajiCopula(copula2 = gumbelCopula(),
                             shapes = fixParam(c(NA_real_, 1),
                                               c(FALSE, TRUE))),
             x = Xdan, start = c(1.1, 0.5), optim.method = "Nelder-Mead")
)

k <- 50
set.seed(7)
withTime(xvCopula(joeCopula(), x = Xdan, k = k))

set.seed(13)
withTime(xvCopula(gumbelCopula(), x = Xdan, k = k))

set.seed(14)
withTime(
    xvCopula(khoudrajiCopula(copula2 = gumbelCopula(),
                             shapes = fixParam(c(NA_real_, 1),
                                               c(FALSE, TRUE))),
             x = Xdan, k = k, start = c(1.1, 0.5),
             optim.method = "Nelder-Mead")
)


### Cross-validation for the rdj data set

withTime(xvCopula(normalCopula(dim = 3, dispstr = "un"), x = Xrdj))

withTime(xvCopula(tCopula(dim = 3, dispstr = "un", df = 10, df.fixed = TRUE),
                  x = Xrdj))

set.seed(22)
withTime(xvCopula(normalCopula(dim = 3, dispstr = "un"), x = Xrdj, k = k))

set.seed(4)
withTime(xvCopula(tCopula(dim = 3, dispstr = "un", df = 10, df.fixed = TRUE),
                  x = Xrdj, k = k))

set.seed(1980)
withTime(xvCopula(tCopula(dim = 3, dispstr = "un"), x = Xrdj, k = k))

summary(fitCopula(tCopula(dim = 3, dispstr = "un"), data = pobs(Xrdj)))

