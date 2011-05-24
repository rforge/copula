\name{K}
\alias{K}
\title{Kendall Distribution Function of a Copula}
\description{
  Compute the Kendall distribution function of an Archimedean copula,
  defined as
  \deqn{K(t) = P(C(U_1,U_2,\dots,U_d) \le t),}{K(t) = P(C(U1,U2,\dots,Ud) <= t),}
  where \eqn{t \in [0,1]}{t in [0,1]}, and the \eqn{d}-dimensional
  \eqn{(U_1,U_2,\dots,U_d)}{(U1,U2,\dots,Ud)} is distributed according
  to the copula \eqn{C}.

  Note that this \dQuote{distribution function of the probability
  integral transform} is known to be the uniform \eqn{K(t) = t}, for
  \eqn{d = 1}, but is non-trivial for \eqn{d \ge 2}{d >= 2}.
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
t <- seq(0,1, length.out = 255)
Kt <- K(t, cop=cop@copula, d=d) # exact
Kt.MC <- K(t, cop=cop@copula, d=d, n.MC=1000) # via Monte Carlo

## build sample from K
set.seed(1)
U <- rnacopula(n, cop)
W <- pnacopula(cop, u=U)

## plot empirical distribution function based on W
## and the corresponding theoretical Kendall distribution function
## (exact and via Monte Carlo)
plot(ecdf(W), col="blue", xlim=c(0,1), verticals=TRUE,
     main = expression("Empirical"~ F[n]( C(U) ) ~
                       "and its Kendall distribution"~K(t)),
     do.points=FALSE, asp=1)
abline(0,1, lty=2); abline(h=0:1, v=0:1, lty=3, col="gray")
lines(t, Kt.MC, col="red")
lines(t, Kt, col="black")
legend(.2,.6, expression(F[n],K(t), K[MC](t)),
       col=c("blue","red","black"), lty=1, bty="n",
       xjust = 0.25, yjust = 0.5)
}
\keyword{distribution}
