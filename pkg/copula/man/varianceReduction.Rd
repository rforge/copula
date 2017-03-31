\name{varianceReduction}
\alias{rAntitheticVariates}
\alias{rLatinHypercube}
\title{Variance-Reduction Methods}
\description{
  Computing antithetic variates or Latin hypercube samples.
}
\usage{
rAntitheticVariates(u)
rLatinHypercube(u, ...)
}
\arguments{
  \item{u}{\eqn{n\times d}{n x d}-\code{\link{matrix}} (or \eqn{d}-\code{\link{vector}}) of random
    variates in the unit hypercube.}
  \item{\dots}{additional arguments passed to the underlying \code{\link{rank}()}.}
}
\value{
  \describe{
    \item{\code{rAntitheticVariates()}}{\code{\link{array}} of dimension
    \eqn{n\times d\times 2}{n x d x 2}, where the last component being 1
      contains the original sample \code{u} and the last component being
      2 contains the sample \code{1-u}.}
    \item{\code{rLatinHypercube()}}{\code{\link{matrix}} of the same
      dimensions as \code{u}.}
  }
}
\details{
  \code{rAntitheticVariates()} takes any copula sample \eqn{u}, builds
  \eqn{1-u}, and returns the two matrices in the form of an array; this
  can be used for the variance-reduction method of (componentwise)
  antithetic variates.

  \code{rLatinHypercube()} takes any copula sample, componentwise
  randomizes its ranks minus 1 and then divides by the sample size in
  order to obtain a Latin hypercubed sample.
}
%\author{Marius Hofert}
%% \seealso{
%% }
\references{
  Cambou, M., Hofert, M. and Lemieux, C. (2016).
  Quasi-random numbers for copula models.
  \emph{Statistics and Computing}, 1--23.

  Packham, N. and Schmidt, W. M. (2010).
  Latin hypercube sampling with dependence and applications in finance.
  \emph{Journal of Computational Finance} \bold{13}(3), 81--111.
}
\examples{
### 1 Basic plots ##############################################################

## Generate data from a Gumbel copula
cop <- gumbelCopula(iTau(gumbelCopula(), tau = 0.5))
n <- 1e4
set.seed(271)
U <- rCopula(n, copula = cop)

## Transform the sample to a Latin Hypercube sample
U.LH <- rLatinHypercube(U)

## Plot
## Note: The 'variance-reducing property' is barely visible, but that's okay
layout(rbind(1:2))
plot(U,    xlab = expression(U[1]), ylab = expression(U[2]), pch = ".")
plot(U.LH, xlab = expression(U[1]), ylab = expression(U[2]), pch = ".")
layout(1) # reset layout

## Transform the sample to an Antithetic variate sample
U.AV <- rAntitheticVariates(U)

## Plot original sample and its corresponding (componentwise) antithetic variates
layout(rbind(1:2))
plot(U.AV[,,1], xlab = expression(U[1]), ylab = expression(U[2]), pch = ".")
plot(U.AV[,,2], xlab = expression(U[1]), ylab = expression(U[2]), pch = ".")
layout(1) # reset layout


### 2 Small variance-reduction study for exceedance probabilities ##############

## Auxiliary function for approximately computing P(U_1 > u_1,..., U_d > u_d) by
## Monte Carlo simulation based on pseudo-sampling and Latin hypercube sampling
survival_prob <- function(n, copula, u = rep(0.99, dim(copula)))
{
    d <- length(u)
    stopifnot(n >= 1, inherits(copula, "Copula"), 0 < u, u < 1, d == dim(copula))
    ## Pseudo-sampling
    U  <- rCopula(n, copula = copula)
    prob <- mean(rowSums(U > rep(u, each = n)) == d)
    ## Latin hypercube sampling
    U. <- rLatinHypercube(U)
    prob. <- mean(rowSums(U. > rep(u, each = n)) == d)
    ## Return
    c(PRNG = prob, LHS = prob.)
}

## Simulate the probabilities of falling in (u_1,1] x ... x (u_d,1].
B <- 500 # number of replications
n <- 1000 # sample size
d <- 2 # dimension
u <- rep(0.99, d) # lower-left endpoint of the considered cube
nu <- 3 # degrees of freedom
th <- iTau(tCopula(df = nu), tau = 0.5) # correlation parameter
cop <- tCopula(param = th, dim = d, df = nu) # t copula
set.seed(271) # for reproducibility
system.time(res <- lapply(1:B, function(b) survival_prob(n, copula = cop, u = u)))

## Grab out the values
PRNG <- sapply(res, function(x) x[["PRNG"]])
LHS  <- sapply(res, function(x) x[["LHS"]])

## Compute the variance-reduction factor and % improvement
(vrf <- var(PRNG) / var(LHS)) # estimated variance reduction factor
(pim <- (var(PRNG) - var(LHS)) / var(PRNG) * 100) # % improvement

## Boxplot
boxplot(list(PRNG = PRNG, LHS = LHS),
        sub = substitute("Variance-reduction factor"~v~
                         "corresponding to an improvement of"~p*"\%",
                         list(v = round(vrf, 2), p = round(pim, 2))),
        main = substitute("Simulated exceedance probabilities"~
                          P(bold(U) > bold(u))~"for a"~t[nu.]~"copula",
                          list(nu. = nu)))
lab <- substitute("Based on B ="~B.~"replications with n ="~n.~"and d ="~d.,
                  list(B. = B, n. = n, d. = d))
mtext(lab, side = 4, line = 1, adj = 0)
}
\keyword{distribution}
\keyword{methods}

