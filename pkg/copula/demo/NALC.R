### Simulating dependent multivariate Levy processes based on positive (nested)
### Archimedean Levy copulas (here: Clayton)

require(gsl) # for exponential integral
require(copula)

doPDF <- Sys.getenv("USER")=="mhofert"


## Auxiliary functions #########################################################

## Tail integral of a variance gamma Levy process
## \bar{\nu}(x) = \int_x^\infty f(z) dz for the Levy density
## f(z) = (c/x)*exp(-lambda*x) for x>0, c=1/kappa and
## lambda=(sqrt(theta^2+2*sigma^2/kappa)-theta)/sigma^2
nu_bar_vargamma <- function(x, theta, kappa, sigma) {
    lambda <- (sqrt(theta^2+2*sigma^2/kappa)-theta)/sigma^2
    -expint_Ei(-lambda*x, give=FALSE)/kappa
}

## Test
if(FALSE) {
    theta <- -0.2
    kappa <- 0.05
    sigma <- 0.3
    x <- seq(.Machine$double.xmin, 29, by=0.01)
    y <- nu_bar_vargamma(x, theta=theta, kappa=kappa, sigma=sigma) # monotonically decreasing
    range(y) # => max. value ~ 14093
    plot(x, y, type="l", log="y")
    mtext(substitute("y-range: ("*l*","~u*")", list(l=min(y), u=max(y))), side=3, line=1)
}

## Inverse of the tail integral of a variance gamma Levy process
## \bar{\nu}(x) = \int_x^\infty f(z) dz for the Levy density
## f(z) = (c/x)*exp(-lambda*x) for x>0, c=1/kappa and
## lambda=(sqrt(theta^2+2*sigma^2/kappa)-theta)/sigma^2
nu_bar_inv_vargamma <- function(Gamma, theta, kappa, sigma, ...)
{
    max.val <- nu_bar_vargamma(.Machine$double.xmin, theta=theta, kappa=kappa, sigma=sigma)
    res <- numeric(length(Gamma))
    large <- Gamma >= max.val
    res[large] <- 0 # de facto indistinguishable from 0 anyways
    if(any(!large)) {
        lambda <- (sqrt(theta^2+2*sigma^2/kappa)-theta)/sigma^2
        nu_bar_vargamma_minus <- function(x, z) -expint_Ei(-lambda*x, give=FALSE)/kappa - z
        res[!large] <- vapply(Gamma[!large], function(Gamma.)
            uniroot(nu_bar_vargamma_minus, z=Gamma.,
                    interval=c(.Machine$double.xmin, 29), ...)$root, NA_real_)
    }
    res
}

## V_0 for Clayton Levy copulas
## tau = truncation point
V0_Clayton_Levy <- function(theta, tau)
{
    Z <- numeric(0) # generate arrival times of a Poi(1) process
    repeat {
        E <- rexp(1)
        if(sum(Z) + E > tau) break
        Z <- c(Z, E)
    }
    W <- cumsum(Z) # compute jump times of a Poi(1) process
    # note: length(W) ~ tau
    (W/theta * gamma(1/theta))^theta # compute V_0 = F_0^{-1}(W)
}

## V_{01} for nested Clayton Levy copulas
## Note: V_{01,k} | V_{0,k} ~ LS^{-1}[\bar{\psi}_{01}(.; V_{0,k})] with
##       \bar{\psi}_{01}(t; V_{0,k}) = \exp(-V_{0,k} t^{\theta_0/\theta_1})
##       = copGumbel@V01() (not copClayton@V01()!)
V01_nested_Clayton_Levy <- function(V0, theta0, theta1)
    copGumbel@V01(V0, theta0=theta0, theta1=theta1)

## \bar{\psi} for Clayton Levy copulas
psi_bar_Clayton <- function(t, theta) t^(-1/theta)

## Transforming Gamma with variance-gamma Levy margins
hom_Levy_process <- function(Gamma, theta, kappa, sigma)
{
    U <- runif(nrow(Gamma)) # jump times
    ord <- order(U) # determine the order of the U's
    jump_time <- U[ord] # (sorted) jump times
    jump_size <- apply(Gamma, 2, function(y)
        nu_bar_inv_vargamma(y, theta=theta, kappa=kappa, sigma=sigma)) # (unsorted) jump sizes (apply inverses of marginal tail integrals)
    value <- apply(jump_size, 2, function(x) cumsum(x[ord])) # sort jump sizes according to U's and add them up => (L_t) at jump times
    list(jump_time=jump_time, value=value)
}

## Plot a multivariate Levy process
plot_Levy_process <- function(jump_times, L, ...)
{
    stopifnot(is.matrix(L), (d <- ncol(L)) >= 2,
              length(jump_times)==nrow(L))
    palette <- colorRampPalette(c("black", "royalblue3", "darkorange2", "maroon3"), space="Lab")
    cols <- palette(d) # d colors
    plot(jump_times, L[,1], type="l", ylim=range(L),
         xlab="t", ylab=expression(bold(L)[t]), col=cols[1],
         ...)
    for(j in 2:d)
        lines(jump_times, L[,j], col=cols[j])
    legend("bottomright", bty="n", lty=rep(1, d), col=cols,
           legend=as.expression( lapply(1:d, function(j) bquote(L[list(t,.(j))]))))
}


### 4d positive Clayton Levy copula ############################################

set.seed(271)
theta <- 4 # theta
d <- 4 # dimension
tau <- 2000 # truncation point

## Sampling Gamma
V0 <- V0_Clayton_Levy(theta, tau=tau) # generate V0
E <- matrix(rexp(length(V0)*d), ncol=d) # Exp(1)
Gamma <- psi_bar_Clayton(E/rep(V0, d), theta=theta) # Gamma

## Plot
if(FALSE) {
    plot(Gamma[,1], type="l", ylim=range(Gamma), xlab="k", ylab="")
    lines(Gamma[,2], col="royalblue3")
    lines(Gamma[,3], col="darkorange2")
    lines(Gamma[,4], col="maroon3")
    legend("topleft", bty="n", lty=rep(1, 4),
       col=c("black", "royalblue3", "darkorange2", "maroon3"),
       legend=c(expression(Gamma[list(k,1)]), expression(Gamma[list(k,2)]),
                expression(Gamma[list(k,3)]), expression(Gamma[list(k,4)])))
}

## Sampling (L_t)
L <- hom_Levy_process(Gamma, theta=-0.2, kappa=0.05, sigma=0.3)

## Plot
if(doPDF) pdf(file=(file <- "fig_L_with_positive_Clayton_Levy_copula.pdf"),
              width=7, height=7)
par(pty="s")
plot_Levy_process(L$jump_time, L$value,
                  main="Levy process with positive Clayton Levy copula")
if(doPDF) dev.off.pdf(file)


### 4d positive nested Clayton Levy copula #####################################

set.seed(271)
theta <- c(0.7, 3, 2) # theta_0, theta_1, theta_2
d <- 4 # dimension
tau <- 2000 # truncation point

## Sampling Gamma
V0 <- V0_Clayton_Levy(theta[1], tau=tau) # generate V0
V01 <- V01_nested_Clayton_Levy(V0, theta0=theta[1], theta1=theta[2]) # compute V_{01}
V02 <- V01_nested_Clayton_Levy(V0, theta0=theta[1], theta1=theta[3]) # compute V_{02}
E <- matrix(rexp(length(V0)*d), ncol=d) # Exp(1)
Gamma <- cbind(psi_bar_Clayton(E[,1:2]/V01, theta=theta[2]),
               psi_bar_Clayton(E[,3:4]/V02, theta=theta[3])) # Gamma

## Plot
if(FALSE) {
    plot(Gamma[,1], type="l", ylim=range(Gamma), xlab="k", ylab="")
    lines(Gamma[,2], col="royalblue3")
    lines(Gamma[,3], col="darkorange2")
    lines(Gamma[,4], col="maroon3")
    legend("topleft", bty="n", lty=rep(1, 4),
       col=c("black", "royalblue3", "darkorange2", "maroon3"),
       legend=c(expression(Gamma[list(k,1)]), expression(Gamma[list(k,2)]),
                expression(Gamma[list(k,3)]), expression(Gamma[list(k,4)])))
}

## Sampling (L_t)
L <- hom_Levy_process(Gamma, theta=-0.2, kappa=0.05, sigma=0.3)

## Plot
if(doPDF) pdf(file=(file <- "fig_L_with_positive_nested_Clayton_Levy_copula.pdf"),
              width=7, height=7)
par(pty="s")
plot_Levy_process(L$jump_time, L$value,
                  main="Levy process with positive nested Clayton Levy copula")
if(doPDF) dev.off.pdf(file)
