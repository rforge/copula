\name{emle}
\title{Maximum Likelihood Estimators for (Nested) Archimedean Copulas}
\alias{emle}
\alias{.emle}
\description{
  Compute (simulated) maximum likelihood estimators for (nested)
  Archimedean copulas.
}
\usage{
emle(u, cop, n.MC=0, optimizer="optimize", method,
     interval=paraOptInterval(u, cop@copula@name),
     start=list(theta=mean(interval)),
     \dots)
.emle(u, cop, n.MC=0,
      interval=paraOptInterval(u, cop@copula@name), \dots)
}
\arguments{
  \item{u}{\eqn{n\times d}{n x d}-matrix of (pseudo-)observations (each value in \eqn{[0,1]}) from the copula, where \eqn{n} denotes the sample size and \eqn{d} the dimension.}
  \item{cop}{\code{\linkS4class{outer_nacopula}} to be estimated (currently only Archimedean copulas are provided).}
  \item{n.MC}{\code{\link{integer}}, if positive, \emph{simulated} maximum likelihood
    estimation (SMLE) is used with sample size equal to
    \code{n.MC}; otherwise (\code{n.MC=0}), MLE. In SMLE, the \eqn{d}th generator derivative and thus the copula density is evaluated via
    (Monte Carlo) simulation, whereas MLE uses the explicit formulas for the generator derivatives as found by Hofert et al. (2011b).}
  \item{optimizer}{a string or \code{NULL}, indicating the optimizer to
    be used, where \code{NULL} means to use \code{\link{optim}} via the
    standard \R function \code{\link[stats4]{mle}()} from package \pkg{stats4},
    whereas the default, \code{"optimize"} uses \code{\link{optimize}} via
    the \R function \code{\link[bbmle]{mle2}()} from package \pkg{bbmle}.}
  \item{method}{only when \code{optimizer} is \code{NULL}, the method to
    be used for \code{\link{optim}}.}
  \item{interval}{nbivariate vector denoting the interval where
	    optimization takes place. The default is computed as described in Hofert et al. (2011a).}
  \item{start}{\code{\link{list}} of initial values, passed through.}
  \item{\dots}{additional parameters passed to \code{\link{optimize}}.}
}
\details{
   Exact formulas for the generator derivatives were derived in Hofert et al. (2011b). Based on these formulas one can compute the (log-)densities of the Archimedean copulas. Note that for some densities, the formulas are numerically highly non-trivial to compute and considerable efforts were put in to make the computations numerically feasible even in large dimensions (see the source code of the Gumbel copula, for example). Both MLE and SMLE showed good performance in the simulation study conducted by Hofert et al. (2011a) including the challenging one hundred dimensional case. Alternative estimators (see also \code{\link{enacopula}}) often used because of their numerical feasibility, might break down in much smaller dimensions.  
}
\value{
  \describe{
    \item{emle}{
      an \R object of class \code{"\link[bbmle:mle2-class]{mle2}"} (and thus useful for obtaining confidence intervals) with the
      (simulated) maximum likelihood estimator.}
    \item{.emle}{\code{\link{list}} as returned by
      \code{\link{optimize}()} including the maximum likelihood estimator (does not confidence intervals but is typically faster).}
  }
}
\author{Martin Maechler, Marius Hofert.}
\references{
 	Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2011a),
   Estimators for Archimedean copulas in high dimensions: A comparison,
   to be submitted.

   Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2011b),
   Likelihood inference for Archimedean copulas,
   to be submitted.
}
\seealso{
  \code{\link[bbmle]{mle2}} from package \pkg{bbmle} and
  \code{\link[stats4]{mle}} from \pkg{stats4} on which \code{mle2} is
  modeled. \code{\link{enacopula}} (wrapper for different estimators).
}
\examples{
tau <- 0.25
(theta <- copGumbel@tauInv(tau)) # 4/3
d <-  20
(cop <- onacopulaL("Gumbel", list(theta,1:d)))

set.seed(1)
n <- 200
U <- rnacopula(n,cop)

## Estimation
system.time(efm <- emle(U, cop))
summary(efm)

## Profile likelihood plot
pfm <- profile(efm)
(ci <- confint(pfm, level=0.95))
stopifnot(ci[1] <= theta, theta <= ci[2])
plot(pfm)               # |z| against theta (|z| = square root of the deviance)
plot(pfm, absVal=FALSE, #  z  against theta 
      show.points=TRUE) # showing how it's interpolated
## and show the true theta:
abline(v=theta, col="lightgray", lwd=2, lty=2)
axis(1, pos = 0, at=theta, label=expression(theta[0]))

## Plot of the log-Likelihood 
logL <- function(x) -efm@minuslogl(x) # -sum(copGumbel@dacopula(U, theta=x, log=TRUE))
logL. <- Vectorize(logL)
I <- paraOptInterval(U, "Gumbel")
curve(logL., from=I[1], to=I[2], xlab="theta", ylab="log-likelihood")
abline(v=theta, col="lightgray", lwd=2, lty=2)
axis(1, pos = 0, at=theta, label=expression(theta[0]))
abline(v=ci[1], col="magenta", lwd=2, lty=2)
abline(v=ci[2], col="magenta", lwd=2, lty=2)
abline(v=efm@coef, col="lightgray", lwd=2, lty=2)
axis(1, pos = 0, at=efm@coef, label=expression(hat(theta)[n]))
}
\keyword{models}
