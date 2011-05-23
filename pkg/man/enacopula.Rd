\name{enacopula}
\alias{enacopula}
\title{Estimation Procedures for Nested Archimedean Copulas}
\description{
  A set of ten different estimators, currently for one-parameter
  Archimedean copulas, of possibly quite high dimensions.
}
\usage{
enacopula(u, cop,
          method = c("mle", "smle", "dmle",
            "mde.normal.CvM", "mde.normal.KS", "mde.log.CvM", "mde.log.KS",
             "tau.tau.mean", "tau.theta.mean", "beta"),
          n.MC = if (method == "smle") 10000 else 0,
          interval = paraOptInterval(u, cop@copula@name),
          xargs = list(), \dots)
}
\arguments{
  \item{u}{\eqn{n\times d}{n x d}-matrix of (pseudo-)observations (each value in \eqn{[0,1]}) from the copula, where \eqn{n} denotes the sample size and \eqn{d} the dimension. Consider applying the function \code{\link{pobs}} first in order to obtain \code{u}.}
  \item{cop}{\code{\linkS4class{outer_nacopula}} to be estimated (currently only Archimedean copulas are provided).}
  \item{method}{a \code{\link{character}} string specifying the estimation method to be used, which has to be one (or a unique abbreviation) of
    \describe{
      \item{\code{"mle"}}{maximum likelihood estimator (MLE) computed via \code{\link{emle.}}.}
      \item{\code{"smle"}}{simulated maximum likelihood estimator (SMLE) computed with the function \code{\link{emle.}}, where \code{n.MC} gives the Monte Carlo sample size.}
      \item{\code{"dmle"}}{MLE based on the diagonal (DMLE), see \code{\link{edmle}}.}
      \item{\code{"mde.normal.CvM"}}{minimum distance estimator based
	on the chisq distribution and Cramer-von Mises distance, see \code{\link{emde}}.}
      \item{\code{"mde.normal.KS"}}{minimum distance estimation based on
	the chisq distribution and Kolmogorov-Smirnov distance, see \code{\link{emde}}.}
      \item{\code{"mde.log.CvM"}}{minimum distance estimation based on
	the Erlang distribution and Cramer-von Mises distance, see \code{\link{emde}}.}
      \item{\code{"mde.log.KS"}}{minimum distance estimation based on
	the Erlang distribution and Kolmogorov-Smirnov distance, see \code{\link{emde}}.}
      \item{\code{"tau.tau.mean"}}{averaged pairwise Kendall's tau estimator}
      \item{\code{"tau.theta.mean"}}{average of pairwise Kendall's tau estimators}
      \item{\code{"beta"}}{multivariate Blomqvist's beta estimator}
    }
  }
  \item{n.MC}{\code{\link{integer}}, if positive, simulated maximum likelihood estimation (SMLE) is used with sample size equal to \code{n.MC}.}
  \item{interval}{bivariate vector denoting the interval where
	    optimization takes place. The default is computed as described in Hofert et al. (2011a). Used for all methods except \code{"tau.tau.mean"} and \code{"tau.theta.mean"}.}
  \item{xargs}{list of additional arguments for the chosen estimation method.}
  \item{\dots}{additional arguments passed to \code{\link{optimize}}.}
}
\details{
 	\code{\link{enacopula}} serves as a wrapper for the different implemented estimators and provides a uniform framework to utilize them. For more information, see the single estimators as given in the section \dQuote{See Also}. 
}
\value{
  the estimated parameter, \eqn{\hat{\theta}}{hat(theta)}, i.e., currently a
  number as only one-parameter Archimedean copulas are considered.
}
\author{Marius Hofert and Martin Maechler}
\references{
  	Hofert, M., Maechler, M., and McNeil, A. J. (2011a),
  \emph{Estimators for Archimedean copulas in high dimensions: A comparison},
   to be submitted.
}
\seealso{
  \code{\link{emle}} which returns an object of
  \code{"\linkS4class{mle}"} providing useful methods not
  available for other estimators. \code{\link{edmle}} for the diagonal maximum likelihood estimator. \code{\link{emde}} for the minimum distance estimators.
  \code{\link{etau}} for the estimators based on Kendall's tau.
  \code{\link{ebeta}} for the estimator based on Blomqvist's beta.
}
\examples{
n <- 200
d <- 20
tau <- 0.25
(theta <- copGumbel@tauInv(tau)) # 4/3
(cop <- onacopulaL("Gumbel", list(theta,1:d)))

set.seed(1)
U <- rnacopula(n, cop)

(meths <- eval(formals(enacopula)$method)) # "mle", "smle", ... 
fun <- function(meth, u, cop, theta){
	run.time <- system.time(val <- enacopula(u, cop=cop, method=meth))
	list(value=val, error=val-theta, utime.ms=1000*run.time[[1]])
}
(res <- sapply(meths, fun, u=U, cop=cop, theta=theta))
}
\keyword{models}
