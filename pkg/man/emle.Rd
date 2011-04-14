\name{emle}
\title{Maximum Likelihood Estimation for (Nested) Archimedean Copulas}
\alias{emle}
\alias{.emle}
\description{
  Compute (simulated) maximum likelihood estimators for (nested)
  Archimedean copulas.
}
\usage{
emle(u, cop, n.MC=0, optimizer="optimize", method,
     interval=paraOptInterval(u, cop@copula@name), \dots)

.emle(u, cop, n.MC=0,
      interval=paraOptInterval(u, cop@copula@name), \dots)
}
\arguments{
  \item{u}{matrix of realizations following the copula}
  \item{cop}{outer_nacopula to be estimated}
  \item{n.MC}{integer, if positive, \emph{simulated} maximum likelihood
    estimation (\dQuote{SMLE}) is used with sample size equal to
    \code{n.MC}; otherwise, MLE.  In SMLE, the density is evaluated via
    (\dQuote{Monte Carlo}) simulation, whereas MLE uses high order
    derivatives of the copula generator \eqn{\psi}.}
  \item{optimizer}{a string or \code{NULL}, indicating the optimizer to
    be used, where \code{NULL} means to use \code{\link{optim}} via the
    standard \R \code{\link[stats4]{mle}()} from package \pkg{stats4},
    whereas the default, \code{"optimize"} uses \code{\link{optimize}} via
    \R \code{\link[bbmle]{mle2}()} from package \pkg{bbmle}.}
  \item{method}{only when \code{optimizer} is \code{NULL}, the method to
    be used for \code{\link{optim}}.}
  \item{interval}{numeric vector of length two, denoting the interval
    where optimization takes place.}
  \item{\dots}{additional parameters passed to \code{\link{optimize}}.}
}
\author{Martin Maechler and Marius Hofert.
}
\details{
... ...
}
\value{
  \describe{
    \item{emle():}{
      an \R object of class \code{"\link[bbmle:mle2-class]{mle2}"} with the
      (simulated) maximum likelihood estimator.}
    \item{.emle():}{the \code{\link{list}} returned from
      \code{\link{optimize}()}.}
  }
}
\references{
 ...
}
\seealso{
  \code{\link[bbmle]{mle2}} from package \pkg{bbmle} and
  \code{\link[stats4]{mle}} from \pkg{stats4} on which \code{mle2} is
  modeled.
}
\examples{
n <- 200
d <-  80
tau <- 0.2
(theta <- copJoe@tauInv(tau))# 1.44381
(cop <- onacopulaL("Joe",list(theta,1:d)))

set.seed(1)
U1 <- rnacopula(n,cop)
efm <- emle(U1, cop)
summary(efm)

pfm <- profile(efm)
(ci <- confint(pfm))
stopifnot(ci[1] <= theta, theta <= ci[2])
plot(pfm)
}
\keyword{models}
