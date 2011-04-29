\name{enacopula}
\alias{enacopula}
\title{Estimation Procedures for Nested Archimedean Copulas}
\description{
  A set of ten different estimators, currently for one-dimensional
  (unnested) Archimedean copulas, of possibly quite high dimensions.
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
  \item{u}{numeric (pseudo data) matrix, of values in \eqn{[0,1]}.  Consider
    using \code{enacopula(\link{pobs}(x), ....)} for your data \code{x}.}
  \item{cop}{\code{\linkS4class{outer_nacopula}} to be estimated.
    Currently, only Archimedean (unnested) copulas can be estimated.}
  \item{method}{the estimation method, a \code{\link{character}}
    string. One (or a unique abbreviation) of
    \describe{
      \item{\code{"mle"}}{MLE}
      \item{\code{"smle"}}{SMLE}
      \item{\code{"dmle"}}{MLE based on the diagonal}
      \item{\code{"mde.normal.CvM"}}{minimum distance estimation based
	on the chisq distribution and CvM distance}
      \item{\code{"mde.normal.KS"}}{minimum distance estimation based on
	the chisq distribution and KS distance}
      \item{\code{"mde.log.CvM"}}{minimum distance estimation based on
	the Erlang distribution and CvM distance}
      \item{\code{"mde.log.KS"}}{minimum distance estimation based on
	the Erlang distribution and KS distance}
      \item{\code{"tau.tau.mean"}}{averaged pairwise Kendall's tau estimator}
      \item{\code{"tau.theta.mean"}}{average of Kendall's tau estimators}
      \item{\code{"beta"}}{multivariate Blomqvist's beta estimator}
    }
  }
  \item{n.MC}{integer; if positive, denotes the sample size for SMLE;
    otherwise MLE}%% is it used at all when method != "smle" ??
  \item{interval}{numeric vector of length two denoting the interval where
    optimization takes place, for "mle", "smle", and "dmle".}
  \item{xargs}{list of additional arguments for the estimation
    procedure.}
  \item{\dots}{possibly additional arguments for
    \code{\link{optimize}}.}
}
% \details{
%  TODO ?
% }
\value{
  the estimated parameter, \eqn{\hat{\theta}}{th^}, i.e. a currently a
    number as only one-dimensional Archimedean families are considered.
}
\author{Marius Hofert and Martin Maechler}
\references{
  Hofert, Maechler, McNeil (2011) .... .... %% FIXME
}
\seealso{
  \code{\link{emle}} which returns a an object of
  \code{"\linkS4class{mle}"} which has useful methods not
  available for other estimators.

  \code{\link{emde}}, for minimum distance estimators, and, e.g.,
  \code{\link{ebeta}} for the beta, tau, and diagonal-ML estimators.
}
\examples{
tau <- 0.4
(theta <- copJoe@tauInv(tau))# 2.219
d <-  20
(cop20 <- onacopulaL("Joe",list(theta,1:d)))

(E.meths <- eval(formals(enacopula)$method))
## "mle"  "smle" ... (all 10)

set.seed(7)
U1 <- rnacopula(256, cop20)

system.time( ## takes a few seconds (8--9 on nb-mm3)
 all.est <- sapply(E.meths, function(MET)
               enacopula(U1, cop20, method = MET))
)
cbind(all.est, err = theta - all.est)# mle is best, even in this one example
}
\keyword{models}
