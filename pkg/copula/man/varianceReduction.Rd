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
}
\keyword{distribution}
\keyword{methods}

