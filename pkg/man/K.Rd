\name{K}
\alias{K}
\title{Kendall distribution function}
\description{
  Kendall distribution function, i.e., the distribution function of the probability integral transformation.
}
\usage{
K(t, cop, d, n.MC=0)
}
\arguments{
  \item{t}{evaluation point(s) (have to be in \eqn{[0,1]}).}
  \item{cop}{acopula with specified parameter.}
  \item{d}{dimension.}
  \item{n.MC}{\code{\link{integer}}, if positive, a Monte Carlo approach is applied with sample size equal to \code{n.MC}; otherwise (\code{n.MC=0}) the exact formula is used based on the generator derivatives as found by Hofert et al. (2011b).}
}
\details{
todo
}
\value{Kendall distribution function at t}
\author{Marius Hofert}
\references{
   Hofert, M., Mächler, M., and McNeil, A. J. (2011b),
   \emph{Likelihood inference for Archimedean copulas},
   to be submitted.
}
\seealso{
  \code{\link{gnacopula}} or \code{\link{emde}} (where \code{\link{K}} is used).
}
\examples{
n <- 200
d <- 20
tau <- 0.5
(theta <- copGumbel@tauInv(tau)) # 2
(cop <- onacopulaL("Gumbel", list(theta,1:d)))

## compute Kendall distribution function
t <- (0:n)/n
K. <- K(t, cop=cop@copula, d=d) # exact
K.MC <- K(t, cop=cop@copula, d=d, n.MC=1000) # via Monte Carlo

## build sample from K
set.seed(1)
U <- rnacopula(n, cop)
W <- pnacopula(cop, u=U)

## plot empirical distribution function based on W
## and the corresponding theoretical Kendall distribution function
## (exact and via Monte Carlo)
plot(ecdf(W), col="blue", xlim=c(0,1), verticals=TRUE, do.points=FALSE)
lines(t, K.MC, col="green")
lines(t, K., col="black")
}
\keyword{distribution}
