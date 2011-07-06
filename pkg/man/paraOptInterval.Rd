\name{paraOptInterval}
\alias{paraOptInterval}
\title{Initial Interval or Value for Archimedean Copula Estimation}
\description{
  Compute an initial interval or initial value for optimization/estimation
  routines (only a heuristic; if this fails, choose your own interval or value).
}
\usage{
paraOptInterval(u, family, h=0.15, tau.range=NULL)
}
\arguments{
  \item{u}{matrix of realizations following the copula family specified by \code{family}.}
  \item{family}{Archimedean family to find an initial interval for.}
  \item{h}{distance measured with the concordance measure Kendall's tau.}
  \item{tau.range}{vector containing lower and upper admissible Kendall's tau.}
}
\details{
   The basic idea is to use the explicit and thus very fast diagonal maximum-likelihood estimator for Gumbel's family to find initial values or intervals (for more sophisticated estimation methods).  Given this estimator \eqn{\hat{\theta}^\mathrm{G}}{hat(theta)^G}, the corresponding Kendall's tau is \eqn{\tau^\mathrm{G}(\hat{\theta}^\mathrm{G})}{tau^G(hat(theta)^G)} where \eqn{\tau^\mathrm{G}(\theta)=(\theta-1)/\theta}{tau^G(theta)=(theta-1)/theta} denotes Kendall's tau for Gumbel.  This provides an estimator of Kendall's tau which is typically much faster to evaluate than (pairwise) Kendall's tau.  Given this estimated \sQuote{amount of concordance}, one can obtain an initial value for the provided family via \deqn{\tau^{-1}(\tau^\mathrm{G}(\hat{\theta}^\mathrm{G})),}{tau^(-1)(tau^G(hat(theta)^G)),} where \eqn{\tau^{-1}}{tau^(-1)} denotes the inverse of Kendall's tau of the family for which the initial value is to be computed.  Introducing a \sQuote{distance to concordance} \eqn{h} (provided by \code{h}) and using lower and upper attainable bounds \eqn{\tau_l}{tau_l} and \eqn{\tau_u}{tau_u} (provided by \code{tau.range}, see the implementation for the default values for the implemented families), an initial interval can be constructed via \deqn{[\tau^{-1}(\max\{\tau^\mathrm{G}(\hat{\theta}^\mathrm{G})-h,\tau_l\}),\tau^{-1}(\min\{\tau^\mathrm{G}(\hat{\theta}^\mathrm{G})+h,\tau_u\})].}{[tau^(-1)(max{tau^G(hat(theta)^G)-h,tau_l}),tau^(-1)(min{tau^G(hat(theta)^G)+h,tau_u})].} This basic construction principle was applied by Hofert et al. (2011a) in a large-scale simulation study for various estimators, see the Section \sQuote{See Also}.
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
	\code{\link{enacopula}}, \code{\link{emle}}, \code{\link{edmle}}, \code{\link{emde}}, and \code{\link{ebeta}} (where \code{paraOptInterval} is applied to find initial intervals).
}
\examples{
## Definition of the function:
paraOptInterval

## Generate some data:
tau <- 0.25
(theta <- copGumbel@tauInv(tau)) # 4/3
d <- 20
(cop <- onacopulaL("Gumbel", list(theta,1:d)))

set.seed(1)
n <- 200
U <- rnacopula(n, cop)

## Initial interval:
paraOptInterval(U, family="Gumbel", h=0.15) # contains theta

## Initial value:
paraOptInterval(U, family="Gumbel", h=0) # 1.3195, quite close
}
