---
title: Elements of Copula Modeling with R
author: Marius Hofert, Ivan Kojadinovic, Martin Mächler, Jun Yan
---

## Code from Chapter 2

Below is the R code from Chapter 2 of the book "Elements of Copula Modeling with
R". The code is also available as an [R script](scripts/02_copula.R). Please
[cite](cite.html) the book or package when using the code; in particular, in
publications.

<!-- Copy script here and indent everything by 4 columns -->

    ## By Marius Hofert, Ivan Kojadinovic, Martin Maechler, Jun Yan

    ## R script for Chapter 2 of Elements of Copula Modeling with R


    ### 2.1 Definition and characterization ########################################

    ### Independence copula

    library(copula)
    d <- 2
    ic <- indepCopula(dim = d)

    set.seed(2008)
    u <- runif(d) # a random point in the unit hypercube
    (Pi <- pCopula(u, copula = ic)) # the value of the independence copula at u

    stopifnot(all.equal(Pi, prod(u))) # check

    wireframe2  (ic, FUN = pCopula, # surface plot of the independence copula
                 col.4 = adjustcolor("black", alpha.f = 0.25))
    contourplot2(ic, FUN = pCopula) # contour plot of the independence copula


    ### C-volumes

    a <- c(1/4, 1/2) # lower left end point
    b <- c(1/3, 1) # upper right end point
    stopifnot(0 <= a, a <= 1, 0 <= b, b <= 1, a <= b) # check
    p <- (b[1] - a[1]) * (b[2] - a[2]) # manual computation
    stopifnot(all.equal(prob(ic, l = a, u = b), p)) # check

    n <- 1000 # sample size
    set.seed(271) # set a seed (for reproducibility)
    U <- rCopula(n, copula = ic) # generate a sample of the independence copula
    plot(U, xlab = quote(U[1]), ylab = quote(U[2]))

    set.seed(271)
    stopifnot(all.equal(U, matrix(runif(n * d), nrow = n)))

    n <- 1e6 # a large sample size to obtain a reasonably accurate approximation
    set.seed(314)
    U <- rCopula(n, copula = ic)
    ## Approximate the Pi-volume by the aforementioned proportion
    p.sim <- mean(a[1] < U[,1] & U[,1] <= b[1] & a[2] < U[,2] & U[,2] <= b[2])
    stopifnot(all.equal(p.sim, p, tol = 1e-2))


    ### Frank copula

    d <- 2 # dimension
    theta <- -10 # copula parameter
    fc <- frankCopula(theta, dim = d) # define a Frank copula

    set.seed(2010)
    n <- 5 # sample size
    u <- matrix(runif(n * d), nrow = n) # random points in [0,1]^d
    pCopula(u, copula = fc) # copula values at u

    dCopula(u, copula = fc) # density values at u

    wireframe2(fc, FUN = pCopula, # wireframe plot (copula)
               draw.4.pCoplines = FALSE)
    wireframe2(fc, FUN = dCopula, delta = 0.001) # wireframe plot (density)
    contourplot2(fc, FUN = pCopula) # contour plot (copula)
    contourplot2(fc, FUN = dCopula) # contour plot (density)

    set.seed(1946)
    n <- 1000
    U   <- rCopula(n, copula = fc)
    U0  <- rCopula(n, copula = setTheta(fc, value = 0))
    U10 <- rCopula(n, copula = setTheta(fc, value = 10))
    plot(U,   xlab = quote(U[1]), ylab = quote(U[2]))
    plot(U0,  xlab = quote(U[1]), ylab = quote(U[2]))
    plot(U10, xlab = quote(U[1]), ylab = quote(U[2]))


    ### Clayton copula

    d <- 3
    cc <- claytonCopula(4, dim = d)# theta = 4

    set.seed(2013)
    n <- 5
    u <- matrix(runif(n * d), nrow = n) # random points in the unit hypercube
    pCopula(u, copula = cc) # copula values at u

    dCopula(u, copula = cc) # density values at u

    set.seed(271)
    n <- 1000
    U <- rCopula(n, copula = cc)
    splom2(U, cex = 0.3, col.mat = "black")


    ### Gumbel-Hougaard copula

    d <- 2
    gc <- gumbelCopula(3, dim = d) # theta = 3

    n <- 1000
    set.seed(1993)
    U <- rCopula(n, copula = gc)
    plot(U, xlab = quote(U[1]), ylab = quote(U[2]))
    wireframe2(gc, dCopula, delta = 0.025) # wireframe plot (density)


    ### 2.2 The Frechet-Hoeffding bounds ###########################################

    ### Frechet-Hoeffding bounds

    n <- 100
    set.seed(1980)
    U <- runif(n)
    plot(cbind(U, 1-U), xlab = quote(U[1]), ylab = quote(U[2]))
    plot(cbind(U, U),   xlab = quote(U[1]), ylab = quote(U[2]))

    n.grid <- 26 # number of grid points
    u <- seq(0, 1, length.out = n.grid) # subdivision points in each dimension
    u12 <- expand.grid("u[1]" = u, "u[2]" = u) # build a grid
    W <- pmax(u12[,1] + u12[,2] - 1, 0) # values of W on grid
    M <- pmin(u12[,1], u12[,2]) # values of M on grid
    val.W <- cbind(u12, "W(u[1],u[2])" = W) # append grid
    val.M <- cbind(u12, "M(u[1],u[2])" = M) # append grid
    wireframe2(val.W)
    wireframe2(val.M)
    contourplot2(val.W, xlim = 0:1, ylim = 0:1)
    contourplot2(val.M, xlim = 0:1, ylim = 0:1)


    ### Marshall-Olkin copulas

    ## A Marshall-Olkin copula
    C <- function(u, alpha)
        pmin(u[,1] * u[,2]^(1 - alpha[2]), u[,1]^(1 - alpha[1]) * u[,2])
    alpha <- c(0.2, 0.8)
    val <- cbind(u12, "C(u[1],u[2])" = C(u12, alpha = alpha)) # append C values
    ## Generate data
    n <- 1000 # sample size for scatter plots
    set.seed(712)
    V <- matrix(runif(n * 3), ncol = 3)
    U <- cbind(pmax(V[,1]^(1/(1 - alpha[1])), V[,3]^(1/alpha[1])),
               pmax(V[,2]^(1/(1 - alpha[2])), V[,3]^(1/alpha[2])))
    ## Plots
    wireframe2(val)
    plot(U, xlab = quote(U[1]), ylab = quote(U[2]))


    ### 2.3 Sklar's Theorem ########################################################

    ### First part of Sklar's Theorem - decomposition

    library(mvtnorm)
    d <- 2 # dimension
    rho <- 0.7 # off-diagonal entries of the correlation matrix P
    P <- matrix(rho, nrow = d, ncol = d) # build the correlation matrix P
    diag(P) <- 1
    set.seed(64)
    u <- runif(d) # generate a random evaluation point
    x <- qnorm(u)
    pmvnorm(upper = x, corr = P) # evaluate the copula C at u

    nc <- normalCopula(rho, dim = d)
    pCopula(u, copula = nc) # value of the copula at u

    nu <- 3 # degrees of freedom
    x. <- qt(u, df = nu)
    pmvt(upper = x., corr = P, df = nu) # evaluate the t copula at u

    try( pmvt(upper = x., corr = P, df = 3.5) )

    tc <- tCopula(rho, dim = d, df = nu)
    pCopula(u, copula = tc) # value of the copula at u


    ### Second part of Sklar's Theorem - composition

    H.obj <- mvdc(claytonCopula(4), c("norm", "exp"),
                  list(list(mean = 1, sd = 2), list(rate = 3)))

    set.seed(1979)
    z <- cbind(rnorm(5, mean = 1, sd = 2), rexp(5, rate = 3)) # evaluation points
    pMvdc(z, H.obj) # values of the df at z

    dMvdc(z, H.obj) # values of the corresponding density at z

    set.seed(1975)
    X <- rMvdc(1000, H.obj)
    plot(X, xlab = quote(X[1]), ylab = quote(X[2]))


    ### 2.4 The invariance principle ###############################################

    ### Sampling from a normal or t copula

    n <- 1000 # sample size
    d <- 2 # dimension
    rho <- 0.7 # off-diagonal entry in the correlation matrix P
    P <- matrix(rho, nrow = d, ncol = d) # build the correlation matrix P
    diag(P) <- 1
    nu <- 3.5 # degrees of freedom
    set.seed(271)
    X <- rmvt(n, sigma = P, df = nu) # n multiv. t observations
    U <- pt(X, df = nu) # n ind. realizations from the corresponding copula

    set.seed(271)
    U. <- rCopula(n, tCopula(rho, dim = d, df = nu))
    stopifnot(all.equal(U, U.)) # test of (numerical) equality

    plot(U., xlab = quote(U[1]), ylab = quote(U[2]))
    plot(U, xlab = quote(U[1]), ylab = quote(U[2]))


    ### From a multivariate t distribution to a t copula to a meta-t model

    Y <- qnorm(U) # transform U (t copula) to normal margins
    ind <- c(A = 725, B = 351, C = 734) # use 'plot(X); identify(X)' to find them
    ## Plot function highlighting A, B, C
    plotABC <- function(X, col = adjustcolor("black", 1/2), pch = 19, ...) {
        cols <- adjustcolor(c("red","blue","magenta"), offset = -c(1,1,1,1.5)/4)
        pTxt <- function(j) {
            lab <- names(ind)[j]
            xy <- X[ind[[lab]], , drop = FALSE]
            for(pc. in pch) points(xy, pch = pc., col = cols[j])
            text(xy, label = lab, adj = c(0.5, -0.6), col = cols[j], font = 2)
        }
        par(pty = "s"); plot(X, col = col, asp = 1,...)
        for(j in seq_along(ind)) pTxt(j)
    }
    ## Scatter plot of observations from the multivariate t distribution
    plotABC(X, xlab = quote(X[1]), ylab = quote(X[2]))
    ## Scatter plot of observations from the corresponding t copula
    plotABC(U, xlab = quote(U[1]), ylab = quote(U[2]))
    ## Scatter plot of observations from the meta-t distribution
    plotABC(Y, xlab = quote(Y[1]), ylab = quote(Y[2]))


    ### Verifying the invariance principle

    rho <- 0.6
    P <- matrix(c(1, rho, rho, 1), ncol = 2) # the correlation matrix
    C <- function(u) pCopula(u, copula = normalCopula(rho)) # normal copula
    Htilde <- function(x)
        apply(cbind(log(x[,1]), -log((1-x[,2])/x[,2])), 1, function(x.)
              pmvnorm(upper = x., corr = P))
    qF1tilde <- function(u) exp(qnorm(u))
    qF2tilde <- function(u) 1/(1+exp(-qnorm(u)))
    Ctilde <- function(u) Htilde(cbind(qF1tilde(u[,1]), qF2tilde(u[,2])))
    set.seed(31)
    u <- matrix(runif(5*2), ncol = 2) # 5 random evaluation points
    stopifnot(all.equal(Ctilde(u), C(u)))

    set.seed(721)
    n <- 1000
    X <- rmvnorm(n, mean = c(0,0), sigma = P) # sample from N(0, P)
    ## Sample the copula of X directly
    U <- pnorm(X)
    ## Transform the sample X componentwise
    TX <- cbind(exp(X[,1]), plogis(X[,2])) # note: plogis(x) = 1/(1+exp(-x))
    ## Apply the marginal dfs to get a sample from the copula of TX
    ## Note: qlogis(p) == logit(p) == log(p/(1-p))
    V <- cbind(pnorm(log(TX[,1])), pnorm(qlogis(TX[,2])))
    stopifnot(all.equal(V, U)) # => the samples of the two copulas are the same


    ### 2.5 Survival copulas and copula symmetries #################################

    ### Survival copulas

    cc <- claytonCopula(2)
    set.seed(271)
    U <- rCopula(1000, copula = cc) # sample from the Clayton copula
    V <- 1-U # sample from the survival Clayton copula
    plot(U, xlab = quote(U[1]), ylab = quote(U[2])) # scatter plot
    plot(V, xlab = quote(V[1]), ylab = quote(V[2])) # for the survival copula

    wireframe2(cc,            FUN = dCopula, delta = 0.025)
    wireframe2(rotCopula(cc), FUN = dCopula, delta = 0.025)


    ### Visually assessing radial symmetry and exchangeability

    contourplot2(tCopula(0.7, df = 3.5), FUN = dCopula)
    contourplot2(gumbelCopula(2),        FUN = dCopula)


    ### 2.6 Measures of association ################################################

    ### 2.6.1 Fallacies related to the correlation coefficient #####################

    ### Counterexample to Fallacies 3 and 4

    ## Evaluate the density of C for h_1(u) = 2*u*(u-1/2)*(u-1),
    ## h_2(u) = theta*u*(1-u) and two different thetas
    n.grid <- 20 # number of grid points in each dimension
    u <- seq(0, 1, length.out = n.grid) # subdivision points in each dimension
    u12 <- expand.grid("u[1]" = u, "u[2]" = u) # build a grid
    dC <- function(u, th) 1 + th * (6 * u[,1] * (u[,1]-1) + 1) * (1 - 2*u[,2])
    wireframe2(cbind(u12, "c(u[1],u[2])" = dC(u12, th = -1)))
    wireframe2(cbind(u12, "c(u[1],u[2])" = dC(u12, th =  1)))


    ### Uncorrelatedness versus independence

    n <- 1000
    set.seed(314)
    Z <- rnorm(n)
    U <- runif(n)
    V <- rep(1, n)
    V[U < 1/2] <- -1 # => V in {-1,1}, each with probability 1/2
    X <- cbind(Z, Z*V) # (X_1,X_2)
    stopifnot(cor.test(X[,1], X[,2])$p.value >= 0.05) # H0:`cor=0' not rejected
    Y <- matrix(rnorm(n * 2), ncol = 2) # independent N(0,1)
    ## Plots
    plot(X, xlab = quote(X[1]), ylab = quote(X[2]))
    plot(Y, xlab = quote(Y[1]), ylab = quote(Y[2]))


    ### Counterexample to Fallacy 5

    ## Function to compute the correlation bounds for LN(0, sigma_.^2) margins
    corBoundLN <- function(s, method = c("max", "min")) {
        ## s = (sigma_1, sigma_2)
        if(!is.matrix(s)) s <- rbind(s)
        method <- match.arg(method)
        if(method == "min") s[,2] <- -s[,2]
        (exp((s[,1]+s[,2])^2/2)-exp((s[,1]^2+s[,2]^2)/2)) /
            sqrt(expm1(s[,1]^2)*exp(s[,1]^2)*expm1(s[,2]^2)*exp(s[,2]^2))
    }
    ## Evaluate correlation bounds on a grid
    n.grid <- 20 # number of grid points in each dimension
    s <- seq(0.01, 5, length.out = n.grid) # subdivision points in each dimension
    s12 <- expand.grid("sigma[1]" = s, "sigma[2]" = s) # build a grid
    ## Plots
    wireframe2(cbind(s12, `underline(Cor)(sigma[1],sigma[2])` =
                              corBoundLN(s12, method = "min")))
    wireframe2(cbind(s12, `bar(Cor)(sigma[1],sigma[2])` = corBoundLN(s12)))


    ### 2.6.2 Rank correlation measures ############################################

    ### rho(), iRho(), tau() and iTau()

    theta <- -0.7
    stopifnot(all.equal(rho(normalCopula(theta)), 6 / pi * asin(theta / 2)))
    stopifnot(all.equal(tau(normalCopula(theta)), 2 / pi * asin(theta)))
    theta <- 2
    stopifnot(all.equal(tau(claytonCopula(theta)), theta / (theta + 2)))
    stopifnot(all.equal(tau(gumbelCopula(theta)), 1 - 1 / theta))

    theta <- (0:8)/16
    stopifnot(all.equal(iRho(normalCopula(), rho = 6/pi * asin(theta/2)), theta))
    stopifnot(all.equal(iTau(normalCopula(), tau = 2/pi * asin(theta)),   theta))
    theta <- 1:20
    stopifnot(all.equal(iTau(claytonCopula(), theta / (theta + 2)), theta))
    stopifnot(all.equal(iTau(gumbelCopula(),  1 - 1 / theta),       theta))

    theta <- 3
    iRho(claytonCopula(), rho = rho(claytonCopula(theta)))


    ### Estimating Spearman's rho and Kendall's tau

    rho <- 0.6 # "true" Spearman's rho
    theta <- iRho(claytonCopula(), rho = rho)
    n <- 1000
    set.seed(974); U <- rCopula(n, copula = claytonCopula(theta))
    rho.def <- cor(apply(U, 2, rank))[1,2]      # Spearman's rho manually
    rho.R   <- cor(U, method = "spearman")[1,2] # Spearman's rho from R
    stopifnot(all.equal(rho.def, rho.R)) # the same
    rho.R  # indeed close to 0.6

    tau <- -0.5 # true Kendall's tau
    theta <- iTau(normalCopula(), tau = tau)
    n <- 1000
    set.seed(974)
    U <- rCopula(n, copula = normalCopula(theta))
    p.n <- 0
    for(i in 1:(n-1)) # number of concordant pairs (obviously inefficient)
        for(j in (i+1):n)
            if(prod(apply(U[c(i,j),], 2, diff)) > 0) p.n <- p.n + 1
    tau.def <- 4 * p.n / (n * (n - 1)) - 1   # Kendall's tau manually (slow!)
    tau.R <- cor(U, method = "kendall")[1,2] # Kendall's tau from R
    stopifnot(all.equal(tau.def, tau.R)) # the same
    tau.R # close to "true" -0.5


    ### Spearman's rho and Kendall's tau under counter- and comonotonicity

    n <- 100
    set.seed(75)
    X <- rnorm(n)
    Y <- -X^3 # perfect negative dependence
    rho.counter <- cor(X, Y, method = "spearman")
    tau.counter <- cor(X, Y, method = "kendall")
    stopifnot(rho.counter == -1, tau.counter == -1)
    Z <- exp(X) # perfect positive dependence
    rho.co <- cor(X, Z, method = "spearman")
    tau.co <- cor(X, Z, method = "kendall")
    stopifnot(rho.co == 1, tau.co == 1)


    ### 2.6.3 Tail dependence coefficients #########################################

    ### Four distributions with N(0,1) margins and a Kendall's tau of 0.7

    ## Kendall's tau and corresponding copula parameters
    tau <- 0.7
    theta.n <- iTau(normalCopula(),  tau = tau)
    theta.t <- iTau(tCopula(df = 3), tau = tau)
    theta.c <- iTau(claytonCopula(), tau = tau)
    theta.g <- iTau(gumbelCopula(),  tau = tau)
    ## Samples from the corresponding 'mvdc' objects
    set.seed(271)
    n <- 10000
    N01.marg <- list(list(mean = 0, sd = 1), list(mean = 0, sd = 1))
    ## For the normal copula
    h.n <- mvdc(normalCopula(theta.n), c("norm", "norm"), N01.marg)
    X.n <- rMvdc(n, mvdc = h.n)
    ## For the t copula
    h.t <- mvdc(tCopula(theta.t, df = 3), c("norm", "norm"), N01.marg)
    X.t <- rMvdc(n, mvdc = h.t)
    ## For the Clayton copula
    h.c <- mvdc(claytonCopula(theta.c), c("norm", "norm"), N01.marg)
    X.c <- rMvdc(n, mvdc = h.c)
    ## For the Gumbel-Hougaard copula
    h.g <- mvdc(gumbelCopula(theta.g), c("norm", "norm"), N01.marg)
    X.g <- rMvdc(n, mvdc = h.g)
    ## Plotting information for scatter plots
    a. <- 0.005
    q <- qnorm(c(a., 1-a.)) # quantiles of order a and 1 - a
    lim <- range(q, X.n, X.t, X.c, X.g)
    lim <- c(floor(lim[1]), ceiling(lim[2]))
    ## Function for producing one scatter plot
    plotCorners <- function(X, qu, lim, smooth = FALSE, ...) {
        if(smooth) plot <- function(...) smoothScatter(nrpoints = nrow(X)/4, ...)
        plot(X, xlim = lim, ylim = lim, xlab = quote(X[1]), ylab = quote(X[2]),
             col = adjustcolor("black", 0.5), ...) # or pch = 16
        abline(h = qu, v = qu, lty = 2, col = adjustcolor("black", 0.6))
        ll <- sum(apply(X <= qu[1], 1, all)) * 100 / n
        ur <- sum(apply(X >= qu[2], 1, all)) * 100 / n
        mtext(sprintf("Lower left: %.2f%%, upper right: %.2f%%", ll, ur),
              cex = 0.9, side = 1, line = -1.5)
    }
    cex <- 0.4
    plotCorners(X.n, qu = q, lim = lim, cex = cex)
    plotCorners(X.t, qu = q, lim = lim, cex = cex)
    plotCorners(X.c, qu = q, lim = lim, cex = cex)
    plotCorners(X.g, qu = q, lim = lim, cex = cex)


    ### Computing the coefficients of tail dependence

    ## Clayton copula
    theta <- 3
    lam.c <- lambda(claytonCopula(theta))
    stopifnot(all.equal(lam.c[["lower"]], 2^(-1/theta)),
              all.equal(lam.c[["upper"]], 0))
    ## Gumbel--Hougaard copula
    lam.g <- lambda(gumbelCopula(theta))
    stopifnot(all.equal(lam.g[["lower"]], 0),
              all.equal(lam.g[["upper"]], 2-2^(1/theta)))
    ## Normal copula
    rho <- 0.7
    nu <- 3
    lam.n <- lambda(normalCopula(rho))
    stopifnot(all.equal(lam.n[["lower"]], 0),
              all.equal(lam.n[["lower"]], lam.n[["upper"]]))
    ## t copula
    lam.t <- lambda(tCopula(rho, df = nu))
    stopifnot(all.equal(lam.t[["lower"]],
                        2*pt(-sqrt((nu+1)*(1-rho)/(1+rho)), df = nu + 1)),
              all.equal(lam.t[["lower"]], lam.t[["upper"]]))


    ### Tail dependence of t copulas

    ## Coefficient of tail dependence as a function of rho
    rho <- seq(-1, 1, by = 0.01)
    nu <- c(3, 4, 8, Inf)
    n.nu <- length(nu)
    lam.rho <- sapply(nu, function(nu.) # (rho, nu) matrix
        sapply(rho, function(rho.) lambda(tCopula(rho., df = nu.))[["lower"]]))
    expr.rho <- as.expression(lapply(1:n.nu, function(j)
        bquote(nu == .(if(nu[j] == Inf) quote(infinity) else nu[j]))))
    matplot(rho, lam.rho, type = "l", lty = 1, lwd = 2, col = 1:n.nu,
            xlab = quote(rho), ylab = quote(lambda))
    legend("topleft", legend = expr.rho, bty = "n", lwd = 2, col = 1:n.nu)
    ## Coefficient of tail dependence as a function of nu
    nu. <- c(seq(3, 12, by = 0.2), Inf)
    rho. <- c(-1, -0.5, 0, 0.5, 1)
    n.rho <- length(rho.)
    lam.nu <- sapply(rho., function(rh) # (nu, rho) matrix
        sapply(nu., function(nu) lambda(tCopula(rh, df = nu))[["lower"]]))
    expr <- as.expression(lapply(1:n.rho, function(j) bquote(rho == .(rho.[j]))))
    matplot(nu., lam.nu, type = "l", lty = 1, lwd = 2, col = 1:n.rho,
            xlab = quote(nu), ylab = quote(lambda))
    legend("right", expr, bty = "n", lwd = 2, col = 1:n.rho)


    ### Effect of rho and nu on P(U_1 > u, U_2 > u) for t copulas

    ## Note: All calculations here are deterministic
    n <- 128
    u <- seq(0.95, to = 0.9999, length.out = n) # levels u for P(U_1> u, U_2> u)
    rho <- c(0.75, 0.5) # correlation parameter rho
    nu <- c(3, 4, 8, Inf) # degrees of freedom
    len <- length(rho) * length(nu)
    tail.prob <- matrix(u, nrow = length(u), ncol = 1 + len) # tail probabilities
    expr <- vector("expression", length = len) # vector of expressions
    ltys <- cols <- numeric(len) # line types and colors
    for(i in seq_along(rho)) { # rho
        for(j in seq_along(nu)) { # degrees of freedom
            k <- length(nu) * (i - 1) + j
            ## Create the copula
            cop <- ellipCopula("t", param = rho[i], df = nu[j])
            ## Evaluate P(U_1 > u, U_2 > u) = P(U_1 <= 1 - u, U_2 <= 1 - u)
            tail.prob[,k+1] <- pCopula(cbind(1 - u, 1 - u), copula = cop)
            ## Create plot information
            expr[k] <- as.expression(
                substitute(group("(",list(rho, nu), ")") ==
                           group("(", list(RR, NN), ")"),
                           list(RR = rho[i],
                                NN = if(is.infinite(nu[j]))
                                         quote(infinity) else nu[j])))
            ltys[k] <- length(rho) - i + 1
            cols[k] <- j
        }
    }
    ## Standardize w.r.t. Gauss case
    tail.prob.fact <- tail.prob # for comparison to Gauss case
    tail.prob.fact[,2:5] <- tail.prob[,2:5] / tail.prob[,5]
    tail.prob.fact[,6:9] <- tail.prob[,6:9] / tail.prob[,9]
    ## Plot tail probabilities
    matplot(tail.prob[,1], tail.prob[,-1], type = "l", lwd = 2, lty = ltys,
            col = cols, xlab = quote(P(U[1]>u, U[2]>u)~~"as a function of u"),
            ylab = "")
    legend("topright", expr, bty = "n", lwd = 2, lty = ltys, col = cols)
    ## Plot standardized tail probabilities
    matplot(tail.prob.fact[,1], tail.prob.fact[,-1], log = "y", type = "l",
            lty = ltys, col = cols, lwd = (wd <- 2*c(1,1,1,1.6,1,1,1,1)),
            xlab = quote(P(U[1]>u, U[2]>u)~~
                         "as a function of u standardized by Gauss case"),
            ylab = "")
    legend("topleft", expr, bty = "n", lwd = wd, lty = ltys, col = cols)


    ### Effect of rho and nu on P(U_1 > 0.99, .., U_d > 0.99) for t copulas

    d <- 2:20 # dimensions
    u <- 0.99 # level u where to evaluate P(U_1 > u, ..., U_d > u)
    tail.pr.d <- matrix(d, nrow = length(d), ncol = 1+len)# tail prob; P[,1] = d
    set.seed(271) # set seed due to MC randomness here
    for(i in seq_along(rho)) { # rho
        for(j in seq_along(nu)) { # degrees of freedom
            k <- length(nu) * (i-1) + j
            for(l in seq_along(d)) { # dimension
                ## Create the copula
                cop <- ellipCopula("t", param = rho[i], dim = d[l], df = nu[j])
                ## Evaluate P(U_1 > u,...,U_d > u) = P(U_1 <= 1-u,...,U_d <= 1-u)
                tail.pr.d[l, k+1] <- pCopula(rep(1-u, d[l]), copula = cop)
            }
        }
    }
    ## Standardize w.r.t. Gauss case
    tail.pr.d.fact <- tail.pr.d # for comparison to Gauss case
    tail.pr.d.fact[,2:5] <- tail.pr.d[,2:5] / tail.pr.d[,5]
    tail.pr.d.fact[,6:9] <- tail.pr.d[,6:9] / tail.pr.d[,9]
    ## Plot tail probabilities
    matplot(tail.pr.d[,1], tail.pr.d[,-1], type = "l", log = "y", yaxt = "n",
            lty = ltys, col = cols, lwd = 2, ylab = "",
            xlab = quote(P(U[1] > 0.99, ..., U[d] > 0.99)~~
                         "as a function of d"))
    sfsmisc::eaxis(2, cex.axis = 0.8); axis(1, at = 2)
    legend("topright",   expr[1:4], bty="n", lty=ltys[1:4], col=cols[1:4], lwd=2)
    legend("bottomleft", expr[5:8], bty="n", lty=ltys[5:8], col=cols[5:8], lwd=2)
    ## Plot standardized tail probabilities
    matplot(tail.pr.d.fact[,1], tail.pr.d.fact[,-1], log = "y", type = "l",
            las = 1, lty = ltys, col = cols,
            lwd = (wd <- 2*c(1,1,1,1.6,1,1,1,1)), ylab = "",
            xlab = quote(P(U[1] > 0.99,..,U[d] > 0.99)~~
                         "as a function of d standardized by Gauss case"))
    legend("topleft", expr, bty = "n", lty = ltys, lwd = wd, col = cols)
    axis(1, at = 2)

    ## Joint exceedance probability under the normal copula
    d <- 5
    rho <- 0.5
    u <- 0.99
    set.seed(271)
    ex.prob.norm <- pCopula(rep(1 - u, d), copula = normalCopula(rho, dim = d))
    1 / (260 * ex.prob.norm) # ~ 51.72 years

    ## Joint exceedance probability under the t copula model with 3 df
    ## 1) Via scaling of the probability obtained from the normal copula
    ## Note that the scaling factor was read off from the previous plot
    1 / (2600 * ex.prob.norm) # ~ 5.17 years

    ## 2) Directly using the t copula
    ex.prob.t3 <- pCopula(rep(1 - u, d), copula = tCopula(rho, dim = d, df = 3))
    1 / (260 * ex.prob.t3) # ~ 5.91 years


    ### 2.7 Rosenblatt transform and conditional sampling ##########################

    ### Evaluation of and sampling from C_{j|1,..,j-1}(.|u_1,..,u_{j-1})

    ## Define the copula
    tau <- 0.5
    nu <- 3.5
    theta <- iTau(tCopula(df = nu), tau = tau)
    tc <- tCopula(theta, df = nu)
    ## Evaluate the df C(.|u_1) for several u_1
    u <- c(0.05, 0.3, 0.7, 0.95)
    u2 <- seq(0, 1, by = 0.01)
    ccop <- sapply(u, function(u.)
        cCopula(cbind(u., u2), copula = tc, indices = 2))
    ## Evaluate the function C(u_2|.) for several u_2
    u1 <- seq(0, 1, by = 0.01)
    ccop. <- sapply(u, function(u.)
        cCopula(cbind(u1, u.), copula = tc, indices = 2))

    matplot(ccop, type = "l", lty = 1, lwd = 2,
            col = (cols <- seq_len(ncol(ccop))), ylab = "",
            xlab = substitute(C["2|1"](u[2]~"|"~u[1])~~"as a function of"~
                              u[2]~"for a"~{C^italic(t)}[list(rho,nu)]~"copula",
                              list(nu = nu)))
    legend("bottomright", bty = "n", lwd = 2, col = cols,
           legend = as.expression(lapply(seq_along(u), function(j)
               substitute(u[1] == u1, list(u1 = u[j])))))

    matplot(ccop., type = "l", lty = 1, lwd = 2,
            col = (cols <- seq_len(ncol(ccop.))), ylab = "",
            xlab = substitute(C["2|1"](u[2]~"|"~u[1])~~"as a function of"~
                              u[1]~"for a"~{C^italic(t)}[list(rho,nu)]~"copula",
                              list(nu = nu)))
    legend("center", bty = "n", lwd = 2, col = cols,
           legend = as.expression(lapply(seq_along(u), function(j)
               substitute(u[2] == u2, list(u2 = u[j])))))

    ## Sample from C_{2|1}(.|u_1)
    n <- 1000
    set.seed(271)
    u2 <- runif(n)
    ## Small u_1
    u1 <- 0.05
    U2 <- cCopula(cbind(u1, u2), copula = tc, indices = 2, inverse = TRUE)
    ## Large u_1
    u1. <- 0.95
    U2. <- cCopula(cbind(u1., u2), copula = tc, indices = 2, inverse = TRUE)
    plot(U2, ylab = substitute(U[2]~"|"~U[1]==u, list(u = u1)))
    plot(U2., ylab = substitute(U[2]~"|"~U[1]==u, list(u = u1.)))


    ### Rosenblatt transform

    ## Sample from a Gumbel-Hougaard copula
    gc <- gumbelCopula(2)
    set.seed(271)
    U <- rCopula(1000, copula = gc)
    ## Apply the transformation R_C with the correct copula
    U. <- cCopula(U, copula = gc)
    ## Apply the transformation R_C with the wrong copula
    gc. <- setTheta(gc, value = 4) # larger theta
    U.. <- cCopula(U, copula = gc.)
    plot(U., xlab = quote(U*"'"[1]), ylab = quote(U*"'"[2]))
    plot(U.., xlab = quote(U*"'"[1]), ylab = quote(U*"'"[2]))


    ### Conditional distribution method and quasi-random copula sampling

    ## Define the Clayton copula to be sampled
    cc <- claytonCopula(2)
    ## Sample from the Clayton copula via CDM (based on pseudo-random numbers)
    set.seed(271)
    U.pseudo <- rCopula(1000, copula = indepCopula())
    U.cc.pseudo <- cCopula(U.pseudo, copula = cc, inverse = TRUE)
    ## Sample from the Clayton copula via CDM (based on quasi-random numbers)
    set.seed(271)
    library(qrng)
    U.quasi <- ghalton(1000, d = 2) # sobol() is typically even faster
    U.cc.quasi <- cCopula(U.quasi, copula = cc, inverse = TRUE)
    plot(U.pseudo, xlab = quote(U*"'"[1]), ylab = quote(U*"'"[2]))
    plot(U.quasi, xlab = quote(U*"'"[1]), ylab = quote(U*"'"[2]))
    plot(U.cc.pseudo, xlab = quote(U[1]), ylab = quote(U[2]))
    plot(U.cc.quasi, xlab = quote(U[1]), ylab = quote(U[2]))


    ### Variance reduction

    ## Auxiliary function for approximately computing P(U_1 > u_1,..., U_d > u_d)
    ## by Monte Carlo simulation based on pseudo-random numbers, Latin hypercube
    ## sampling and quasi-random numbers.
    sProb <- function(n, copula, u)
    {
        d <- length(u)
        stopifnot(n >= 1, inherits(copula, "Copula"), 0 < u, u < 1,
                  d == dim(copula))
        umat <- rep(u, each = n)
        ## Pseudo-random numbers
        U <- rCopula(n, copula = copula)
        PRNG <- mean(rowSums(U > umat) == d)
        ## Latin hypercube sampling (based on the recycled 'U')
        U. <- rLatinHypercube(U)
        LHS <- mean(rowSums(U. > umat) == d)
        ## Quasi-random numbers
        U.. <- cCopula(sobol(n, d = d, randomize = TRUE), copula = copula,
                       inverse = TRUE)
        QRNG <- mean(rowSums(U.. > umat) == d)
        ## Return
        c(PRNG = PRNG, LHS = LHS, QRNG = QRNG)
    }
    ## Simulate the probabilities of falling in (u_1, 1] x ... x (u_d, 1]
    library(qrng) # for quasi-random numbers
    B <- 500 # number of replications
    n <- 5000 # sample size
    d <- 5 # dimension
    nu <- 3 # degrees of freedom
    th <- iTau(tCopula(df = nu), tau = 0.5) # correlation parameter
    cop <- tCopula(param = th, dim = d, df = nu) # t copula
    u <- rep(0.99, d) # lower-left endpoint of the considered cube
    set.seed(42) # for reproducibility
    true <- prob(cop, l = u, u = rep(1, d)) # true exceedance probability
    res <- lapply(1:B, function(b) sProb(n, copula = cop, u = u))
    ## Grab out the results
    PRNG <- sapply(res, function(x) x[["PRNG"]])
    LHS  <- sapply(res, function(x) x[["LHS"]])
    QRNG <- sapply(res, function(x) x[["QRNG"]])
    ## Compute the variance-reduction factors and % improvements
    vrf  <- var(PRNG) / var(LHS) # variance-reduction factor w.r.t. LHS
    vrf. <- var(PRNG) / var(QRNG) # variance-reduction factor w.r.t. QRNG
    pim  <- (var(PRNG) - var(LHS))  / var(PRNG) * 100 # improvement w.r.t. LHS
    pim. <- (var(PRNG) - var(QRNG)) / var(PRNG) * 100 # improvement w.r.t. QRNG
    ## Boxplot
    boxplot(list(PRNG = PRNG, LHS = LHS, QRNG = QRNG), sub = substitute(
            "Variance-reduction factors and % improvements:"~v1~"("*p1*"%"*"),"~
            v2~"("*p2*"%)",
            list(v1 = round(vrf, 1), v2 = round(vrf., 1), p1 = round(pim),
                 p2 = round(pim.))),
            main = substitute("Simulated exceedance probabilities"~
                              P(bold(U) > bold(u))~"for a"~t[nu.]~"copula",
                              list(nu. = nu)))
    abline(h = true, lty = 3) # true value
    lab <- substitute("B ="~B.~"replications with n ="~n.~"and d ="~d.,
                      list(B. = B, n. = n, d. = d))
    mtext(lab, side = 4, line = 1, adj = 0, las = 0)

