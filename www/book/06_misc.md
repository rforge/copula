---
title: Elements of Copula Modeling with R
author: Marius Hofert, Ivan Kojadinovic, Martin Mächler, Jun Yan
---

## Code from Chapter 6

Below is the R code from Chapter 6 of the book "Elements of Copula
Modeling with R". The code is also available as an [R
script](scripts/06_misc_code.R).  Please [cite](cite.html) the book or
package when using the code; in particular, in publications.

<!-- Copy script here and indent everything by 4 columns -->

    ## By Marius Hofert, Ivan Kojadinovic, Martin Maechler, Jun Yan

    ## R script for Chapter 6 of Elements of Copula Modeling with R


    ### 6.1 Non-stationarity #######################################################

    ### Test of stationarity based on S_n^H

    data(rdj)
    library(xts)
    Xrdj <- xts(rdj[,2:4], order.by = rdj[,1])

    plot.zoo(Xrdj, main = "", xlab = "", mar = c(0, 7, 0, 2.1))

    library(npcp)
    set.seed(981)
    (res <- withTime( cpTestFn(Xrdj, b = NULL) ))
    ## => test result (component value) and run time (component sys.time)

    out <- res$value # the output of the test only
    rdj[which(out$cvm == out$statistic), 1]

    plot(out$cvm, type = "l", xlab = "k", ylab = quote({S^H}[list(n,k)]))

    data(gasoil) # oil and gas prices
    Rgasoil <- apply(log(gasoil[,-1]), 2, diff) # bivariate daily log-returns
    Xgasoil <- xts(Rgasoil, order.by = gasoil[-1, 1]) # corresponding xts object

    plot.zoo(Xgasoil, main = "", xlab = "", mar = c(0, 7, 0, 2.1))

    set.seed(292); withTime( cpTestFn(Xgasoil, b = NULL) )


    ### Test of stationarity based on S_n^C

    set.seed(314)
    withTime( cpTestCn(Xrdj, b = NULL, method = "nonseq") )

    set.seed(137)
    withTime( cpTestCn(Xgasoil, b = NULL, method = "nonseq") )


    ### 6.2 Serial dependence ######################################################

    ### Correlogram and Ljung-Box test of serial independence

    acf(Xgasoil^2, main = "", ci.col = 1)
    text(c(8, 23.4, 6, 25.2), c(1.095, 1.09, 0.26, 0.26),
         labels = c("Oil", "Oil and gas", "Gas and oil", "Gas"))

    Box.test(Xgasoil[,1]^2, lag =  5, type = "Ljung-Box")

    Box.test(Xgasoil[,2]^2, lag =  5, type = "Ljung-Box")

    Box.test(Xgasoil[,1]^2, lag = 20, type = "Ljung-Box")

    Box.test(Xgasoil[,2]^2, lag = 20, type = "Ljung-Box")


    ### Ljung-Box tests can be too liberal

    pvalBox <- function(n) {
        x <- rt(n, df=4)
        c(lag5  = Box.test(x^2, type = "Ljung-Box", lag =  5)$p.value,
          lag20 = Box.test(x^2, type = "Ljung-Box", lag = 20)$p.value)
    }

    set.seed(3298)
    pv <- replicate(10000, pvalBox(500))
    alpha <- c(0.01, 0.05, 0.1) # nominal levels
    ## Empirical levels
    rbind(nom.level = alpha,
          emp.level.lag5  = ecdf(pv["lag5", ])(alpha),
          emp.level.lag20 = ecdf(pv["lag20",])(alpha))


    ### Tests of serial independence based on S_n^s

    set.seed(137)
    sI.d <- withTime( serialIndepTestSim(nrow(Xgasoil), lag.max = 5) )
    sI.d$sys.time # the run time

    serialIndepTest(Xgasoil[,1]^2, d = sI.d$value)

    serialIndepTest(Xgasoil[,2]^2, d = sI.d$value)


    ### 6.3 Filtering ##############################################################

    ### Filtering based on ARMA-GARCH models

    library(rugarch)
    ## Specify ARMA(1,1)-GARCH(1,1) model with Student t innovations
    meanModel <- list(armaOrder = c(1,1)) # ARMA(1,1)
    varModel  <- list(model = "sGARCH", garchOrder = c(1,1)) # GARCH(1,1)
    uspec <- ugarchspec(varModel, mean.model = meanModel,
                        distribution.model = "std") # Student t

    ## Fit marginal ARMA-GARCH models
    fit <- apply(Xgasoil, 2, function(x) ugarchfit(uspec, data = x))
    ## Extract the standardized residuals and build the pseudo-observations
    Z <- sapply(fit, residuals, standardize = TRUE) # standardized residuals
    (nus <- sapply(fit, function(x) x@fit$coef[["shape"]])) # fitted d.o.f.

    U <- pobs(Z) # pseudo-observations from the residuals Z
    plot(U, xlab = expression(U[1]), ylab = expression(U[2])) # pseudo-obs.

    ## Fit a t copula to the innovations
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
    ##    Note: The innovations Z have to be standardized (mean 0, variance 1)
    ##          for ugarchsim()
    Z. <- sapply(1:2, function(j)
        sqrt((nus[j]-2)/nus[j]) * qt(U.[,j], df = nus[j]))
    ## 3) Feed the (cross-sec. dependent) innovations to the marginal ARMA-GARCH
    ##    models and simulate from them
    sim <- lapply(1:2, function(j)
        ugarchsim(fit[[j]], # fitted marginal ARMA-GARCH model
                  n.sim = n.sim, # sample size
                  m.sim = m.sim, # number of trajectories/paths
                  custom.dist = list(name = "sample", # our innovations
                                     distfit = matrix(Z.[,j], ncol = m.sim))))
    ## 4) Extract the simulated (cross-sec. dependent) series X_t and build the
    ##    corresponding (predicted/simulated) oil and gas prices
    X. <- lapply(sim, function(x) fitted(x)) # equal to seriesSim in x@simulation
    S.t <- as.numeric(tail(gasoil, n = 1)[,2:3]) # last available prices S_t
    library(qrmtools) # for log_returns()
    S. <- lapply(1:2, function(j) # predicted prices
        log_returns(X.[[j]], inverse = TRUE, start = rep(S.t[j], m.sim)))
    S.T <- sapply(1:2, function(j) tail(S.[[j]], n = 1)) # pick out prices at T

    library(MASS)
    pred.dens <- kde2d(S.T[,1], S.T[,2], n = 300, lims = c(0, 200, 0, 15))
    image(pred.dens, xlab = "Oil price", ylab = "Gas price",
          col = gray(seq(1, 0, length.out = 100)))


    ### 6.4 Ties ###################################################################

    ### Computing ranks in the presence of ties

    set.seed(1979)
    (U <- runif(8))

    rank(U) # does not depend on the setting of `ties.method'

    b <- 10 # number of bins
    (U.ties <- cut(U, breaks = 0:b/b, labels = 0.5:(b - 0.5)/b)) # a factor

    rank(U.ties)

    rank(U.ties, ties.method = "max")

    rank(U.ties, ties.method = "min")

    set.seed(8)
    rank(U.ties, ties.method = "random")

    set.seed(46)
    rank(U.ties, ties.method = "random")


    ### Effect of ties on multivariate scaled ranks

    n <- 200
    set.seed(4809)
    U <- rCopula(n, claytonCopula(5))
    b <- 10 # number of bins
    U.ties <- apply(U, 2, cut, breaks = 0:b/b, labels =  0.5:(b - 0.5)/ b)

    ties.method <- "max" # can be changed to "min", "average" or "random"
    stopifnot(all.equal(pobs(U.ties, ties.method = ties.method),
                        apply(U.ties, 2, rank, ties.method = ties.method) /
                        (nrow(U.ties) + 1))) # check

    plot(pobs(U), xlab = "", ylab = "")
    plot(pobs(U.ties), xlim = 0:1, ylim = 0:1, xlab = "", ylab = "")
    set.seed(732)
    plot(pobs(U.ties, ties.method = "random"), xlab = "", ylab = "")

    a <- runif(1, min = 0, max = 1/(2 * b))
    set.seed(732)
    V <- pobs(matrix(as.numeric(U.ties) + runif(2 * n, min = -a, max = a),
                     nrow = n, ncol = 2))
    set.seed(732)
    stopifnot(all.equal(V, pobs(U.ties, ties.method = "random"))) # check


    ### Estimation of copula parameters in the presence of ties

    myCopula <- frankCopula # a one-parameter copula
    d <- 2 # dimension
    tau <- 0.25 # strength of dependence
    theta <- iTau(myCopula(), tau = tau) # corresponding copula parameter

    ## The discretizing function; returns a matrix of real numbers
    discretize <- function(U, b)
        matrix(as.numeric(apply(U,2, cut, breaks = 0:b/b,
                                labels = 0.5:(b - 0.5)/ b)), nrow = nrow(U))
    b <- 10 # number of bins

    optim.method <- "BFGS" # the numerical optimization method for MPLE

    m <- 30 # number of times randomization is done for one discretized sample
    ## Function returning the average of estimates over m randomizations
    ## when performing fitting based on random ranks from one discrete sample
    fitCopulaRan <- function(U.ties, method) {
        oneFitRan <- function() {
            V <- pobs(U.ties, ties.method = "random") # break ties at random
            fitCopula(myCopula(dim = d), data = V, method = method,
                      optim.method = optim.method)@estimate # return param. est.
        }
        mean(replicate(m, oneFitRan())) # average of est. over m randomizations
    }

    oneFit <- function(n) {
        U <- rCopula(n, copula = myCopula(theta, dim = d)) # a sample of size n
        ## Rank-based fitting from the continuous sample
        V <- pobs(U) # pseudo-observations (no ties)
        mpl <- fitCopula(myCopula(dim = d), data = V,
                         optim.method = optim.method)@estimate # MPLE
        itau <- fitCopula(myCopula(dim = d), data = V,
                          method = "itau")@estimate # tau inversion
        ## Fitting from the discretized sample based on average ranks
        U.ties <- discretize(U, b = b) # discretized sample
        W <- pobs(U.ties) # corresponding pseudo-observations (average ranks)
        mpl.ave <- fitCopula(myCopula(dim = d), data = W,
                             optim.method = optim.method)@estimate # MPLE
        itau.ave <- fitCopula(myCopula(dim = d), data = W,
                              method = "itau")@estimate # tau inversion
        ## Fitting from the discretized sample based on random ranks
        mpl.ran <-  fitCopulaRan(U.ties, method = "mpl") # MPLE
        itau.ran <- fitCopulaRan(U.ties, method = "itau") # tau inversion
        ## The different estimates of theta
        theta.n <- c(mpl, mpl.ave, mpl.ran, itau, itau.ave, itau.ran)
        names(theta.n) <- c("mpl", "mpl.ave", "mpl.ran",
                            "itau", "itau.ave", "itau.ran")
        theta.n
    }

    set.seed(2010); est <- withTime( replicate(100, oneFit(n = 100)) )

    sumSims <- function() {
        cat(describeCop(myCopula(), kind = "very short"),
            "with a Kendall's tau of", tau, "\n")
        print(est$sys.time) # run time
        ## Estimates on Kendall's tau scale
        tau.n.s <- apply(est$value, c(1, 2), function(x) tau(myCopula(x)))
        ## Relative biases on Kendall's tau scale
        bias <- (rowMeans(tau.n.s) - tau) / tau
        ## Standard deviations of estimates on Kendall's tau scale
        std <- apply(tau.n.s, 1, sd)
        ## Root mean square errors on Kendall's tau scale
        rmse <- sqrt(rowMeans((tau.n.s - tau)^2))
        round(rbind(bias = bias, std = std, rmse = rmse) * 100, 2)
    }
    sumSims()

    tau <- 0.5; theta <- iTau(myCopula(), tau = tau)
    set.seed(2010); est <- withTime( replicate(100, oneFit(n = 100)) )
    sumSims()

    tau <- 0.75; theta <- iTau(myCopula(), tau = tau)
    set.seed(2010); est <- withTime( replicate(100, oneFit(n = 100)) )
    sumSims()


    ### Tests not adapted to the presence of ties

    tau <- 0.5 # strength of dependence
    theta <- iTau(gumbelCopula(), tau = tau) # corresponding copula parameter
    b <- 10 # number of bins for discretization
    ## Function which applies the four tests to data generated from
    ## a Gumbel-Hougaard copula and then discretized
    pvalTies <- function(n) {
        U <- rCopula(n, copula = gumbelCopula(theta))
        U.ties <- discretize(U, b = b)
        res <- c(exchTest(U.ties, ties = FALSE)$p.value,
                 evTestC(U.ties)$p.value,
                 gofCopula(gumbelCopula(), x = U.ties, estim.method = "itau",
                           ties = FALSE)$p.value, # pb
                 gofCopula(gumbelCopula(), x = U.ties, estim.method = "itau",
                           sim = "mult")$p.value)
        names(res) <- c("exch", "ev", "pb.GH", "mult.GH")
        res
    }

    set.seed(4478)
    pv <- withTime( replicate(100, pvalTies(n = 100)) )

    pv$sys.time

    round(apply(pv$value, 1, max), 4) # largest p-values


    ### Tests of exchangeability, radial symmetry and extreme-value dependence for ties

    ## Function which applies the two tests to data generated from
    ## a Gumbel-Hougaard copula and then discretized
    pvalExchEV <- function(n) {
        U <- rCopula(n, copula = gumbelCopula(theta))
        U.ties <- discretize(U, b = b)
        res <- c(exchTest(U.ties, ties = TRUE)$p.value,
                 evTestK(U.ties, ties = TRUE)$p.value)
        names(res) <- c("exch", "ev")
        res
    }

    set.seed(8848)
    pv <- withTime( replicate(1000, pvalExchEV(n = 100)) )
    pv$sys.time # run time

    alpha <- c(0.01, 0.05, 0.1) # nominal levels
    rbind(nom.level = alpha,
          emp.level.exch = ecdf(pv$value["exch",])(alpha),
          emp.level.ev = ecdf(pv$value["ev",])(alpha))

    myCopula <- frankCopula # a one-parameter radially symmetric copula
    tau <- 0.5 # strength of dependence
    d <- 2 # dimension
    theta <- iTau(myCopula(), tau = tau) # corresponding copula parameter
    b <- 10 # number of bins
    ## Function which applies the tests to data generated from
    ## the copula 'myCopula' and then discretized
    pvalSym <- function(n) {
        U <- rCopula(n, copula = myCopula(theta, dim = d))
        U.ties <- discretize(U, b = b)
        radSymTest(U.ties, ties = TRUE)$p.value
    }

    set.seed(3298)
    pv <- withTime( replicate(1000, pvalSym(n = 100)) )
    pv$sys.time # run time

    rbind(nom.level = alpha, emp.level = ecdf(pv$value)(alpha)) # levels


    ### Parametric bootstrap-based goodness-of-fit test adapted to ties

    myCopula <- gumbelCopula # a one-parameter bivariate copula
    tau <- 0.5 # strength of dependence
    theta <- iTau(myCopula(), tau = tau) # corresponding copula parameter
    b <- 10 # number of bins
    ## Function which applies the parametric goodness-of-fit test to data
    ## generated under the null hypothesis and returns the p-value
    pvalPB <- function(n) {
        U <- rCopula(n, copula = myCopula(theta))
        U.ties <- discretize(U, b = b)
        gofCopula(myCopula(), x = U.ties, estim.method = "itau",
                  ties = TRUE)$p.value
    }

    set.seed(32); pv <- withTime( replicate(1000, pvalPB(n = 100)) )
    pv$sys.time # run time

    alpha <- c(0.01, 0.05, 0.1) # nominal levels
    rbind(nom.level = alpha, emp.level = ecdf(pv$value)(alpha)) # levels

    discretize <- function(U, b)
        matrix(as.numeric(apply(U, 2, cut,
                                breaks = (0:b)^2/b^2,
                                labels =  0.5:(b - 0.5)/ b)), nrow = nrow(U))
    set.seed(32); pv <- withTime( replicate(1000, pvalPB(n = 100)) )

    rbind(nom.level = alpha, emp.level = ecdf(pv$value)(alpha))

    discretize <- function(U, b)
        cbind(U[,1], as.numeric(cut(U[,2], breaks = 0:b/b,
                                    labels = 0.5:(b - 0.5)/b)))
    set.seed(32); pv <- withTime( replicate(1000, pvalPB(n = 100)) )

    rbind(nom.level = alpha, emp.level = ecdf(pv$value)(alpha))

    discretize <- function(U, b) U
    set.seed(32); pv <- withTime( replicate(1000, pvalPB(n = 100)) )

    rbind(nom.level = alpha, emp.level = ecdf(pv$value)(alpha))


    ### Effect of ties on cross-validation

    myCopula <- claytonCopula # a one-parameter bivariate copula
    tau <- 0.5 # strength of dependence
    theta <- iTau(myCopula(), tau = tau) # corresponding copula parameter
    b <- 20 # number of bins
    ## The discretizing function
    discretize <- function(U, b)
        matrix(as.numeric(apply(U, 2, cut, breaks = 0:b/b,
                                labels =  0.5:(b - 0.5)/ b)), nrow = nrow(U))
    optim.method <- "BFGS" # the numerical optimization method for MPLE
    ## 10-fold cross-validation
    ## Returns the label of the copula with the highest xv score
    xv <- function(n) {
        U <- rCopula(n, copula = myCopula(theta))
        U.ties <- discretize(U, b = b)
        score <- c(xvCopula(claytonCopula(), x = U.ties, k = 10,
                            optim.method = optim.method),
                   xvCopula(gumbelCopula(), x = U.ties, k = 10,
                            optim.method = optim.method),
                   xvCopula(frankCopula(), x = U.ties, k = 10,
                            optim.method = optim.method))
        c("C", "GH", "F")[which(score == max(score))]
    }

    set.seed(2885)
    withTime( summary(as.factor(replicate(100, xv(n = 100)))) )

    myCopula <- gumbelCopula
    theta <- iTau(myCopula(), tau = tau)
    set.seed(2885)
    withTime( summary(as.factor(replicate(100, xv(n = 100)))) )

    myCopula <- frankCopula
    theta <- iTau(myCopula(), tau = tau)
    set.seed(2885)
    withTime( summary(as.factor(replicate(100, xv(n = 100)))) )


    ### Analysis of the loss insurance data

    data(loss)
    X <- as.matrix(subset(loss, censored == 0, select = c("loss", "alae")))

    apply(X, 2, function(x) length(unique(x)))

    plot(pobs(X), xlab = "loss", ylab = "alae")
    plot(qnorm(pobs(X)), xlab = "loss", ylab = "alae")

    set.seed(3070)
    withTime(
        round(c(exchTest(X, ties = TRUE)$p.value,
                radSymTest(X, ties = TRUE)$p.value,
                evTestK(X, ties = TRUE)$p.value), 4)
    )

    ## Goodness-of-fit testing
    set.seed(4634)
    withTime(
        round(c(gofCopula(gumbelCopula(), x = X, ties = TRUE,
                          estim.method = "itau")$p.value,
                gofCopula(rotCopula(claytonCopula()), x = X, ties = TRUE,
                          estim.method = "itau")$p.value,
                gofCopula(frankCopula(), x = X, ties = TRUE,
                          estim.method = "itau")$p.value,
                gofCopula(plackettCopula(), x = X, ties = TRUE,
                          estim.method = "itau")$p.value,
                gofCopula(normalCopula(), x = X, ties = TRUE,
                          estim.method = "itau")$p.value), 4)
    )

    ## Model selection
    set.seed(4807)
    k <- 50 # for k-fold cross-validation
    optim.method <- "BFGS" # the numerical optimization method for MPLE
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

    withTime( gofCopula(gumbelCopula(), x = X, ties = FALSE,
                        estim.method = "itau") )

