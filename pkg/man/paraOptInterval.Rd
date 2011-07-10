\name{initOpt}
\alias{initOpt}
\title{Initial Interval or Value for Archimedean Copula Estimation}
\description{
  Compute an initial interval or initial value for optimization/estimation
  routines (only a heuristic; if this fails, choose your own interval or value).
}
\usage{
initOpt(family, tau.range=NULL, value=FALSE, u)
}
\arguments{
  \item{family}{Archimedean family to find an initial interval for.}
  \item{tau.range}{vector containing lower and upper admissible Kendall's tau.}
  \item{value}{logical indicating whether an initial interval (the default) or
    an initial value should be returned.}
  \item{u}{matrix of realizations following the copula family specified by \code{family}.  Note that \code{u} can be omitted if \code{value=FALSE}.}
}
\details{
  The basic idea is to use the explicit and thus very fast diagonal
  maximum-likelihood estimator for Gumbel's family to find initial values or
  intervals (for more sophisticated estimation methods).  Given this estimator
  \eqn{\hat{\theta}^\mathrm{G}}{hat(theta)^G}, the corresponding Kendall's tau
  is \eqn{\tau^\mathrm{G}(\hat{\theta}^\mathrm{G})}{tau^G(hat(theta)^G)} where
  \eqn{\tau^\mathrm{G}(\theta)=(\theta-1)/\theta}{tau^G(theta)=(theta-1)/theta}
  denotes Kendall's tau for Gumbel.  This provides an estimator of Kendall's tau
  which is typically much faster to evaluate than, for example, pairwise
  Kendall's taus.  Given this estimated \sQuote{amount of concordance}, one can
  obtain an initial value for the provided family via
  \deqn{\tau^{-1}(\tau^\mathrm{G}(\hat{\theta}^\mathrm{G})),}{tau^(-1)(tau^G(hat(theta)^G)),} 
  where \eqn{\tau^{-1}}{tau^(-1)} denotes the inverse of Kendall's tau of the
  family for which the initial value is to be computed.  This idea is used to
  compute an initial value if \code{value=TRUE}.  Note that if
  \eqn{\tau^\mathrm{G}(\hat{\theta}^\mathrm{G})}{tau^G(hat(theta)^G)} does not
  lie in the range of Kendall's tau as provided by the bivariate vector
  \code{tau.range}, the point in \code{tau.range} closest to
  \eqn{\tau^\mathrm{G}(\hat{\theta}^\mathrm{G})}{tau^G(hat(theta)^G)} is chosen
  instead.  The default (\code{value=FALSE}) returns a reasonably large initial
  interval, see the default of \code{tau.range} in \code{initOpt} for the chosen
  values (in terms of Kendall's tau).  These default values cover a large range
  of concordance.  If this interval is (still) too small, one can adjust it by
  providing a vector of length two to \code{tau.range}. If it is too large, a
  \sQuote{distance to concordance} can be used to determine parameter values
  such that the corresponding Kendall's taus share a certain distance to the
  initial value.  For more details, see Hofert et al. (2011a).  Finally, let us
  note that \code{u} is not required if \code{value=FALSE}.
}
\value{initial interval which can be used for optimization (e.g., for
  \code{\link{emle}}).
}
\author{Marius Hofert}
\references{
  Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2011a),
  Estimators for Archimedean copulas in high dimensions: A comparison,
  to be submitted.
}
\seealso{
  \code{\link{enacopula}}, \code{\link{emle}}, \code{\link{edmle}}, \code{\link{emde}}, and \code{\link{ebeta}} (where \code{initOpt} is applied to find initial intervals).
}
\examples{
## Definition of the function:
initOpt

## Generate some data:
tau <- 0.25
(theta <- copGumbel@tauInv(tau)) # 4/3
d <- 20
(cop <- onacopulaL("Gumbel", list(theta,1:d)))

set.seed(1)
n <- 200
U <- rnacopula(n, cop)

## Initial interval:
initOpt("Gumbel") # contains theta

## Initial value:
initOpt("Gumbel", value=TRUE, u=U) # 1.3195, quite close
}
