\name{nac-estimators}
\alias{ebeta}
\alias{edmle}
\alias{etau}
\title{Miscellaneous Estimators for Nested Archimedean Copulas}
\description{
  Miscellaneous Estimators for Nested Archimedean Copulas, namely,
  \describe{
    \item{ebeta}{Method-of-moment-like parameter estimation of nested
      Archimedean copulas based on Blomqvist's beta.}
    \item{edmle}{Maximum likelihood estimation based on the diagonal of a nested Archimedean copula.}
    \item{etau}{Compute pairwise estimators for nested Archimedean
      copulas based on Kendall's tau, see \code{\link{cor}(*, method="kendall")}.}
  }
}
\usage{
ebeta(u, cop, interval = paraOptInterval(u, cop@copula@name), ...)
edmle(u, cop, interval = paraOptInterval(u, cop@copula@name), ...)
etau (u, cop, method = c("tau.mean", "theta.mean"), ...)
}
\arguments{
  \item{u}{numeric matrix of realizations following the copula}
  \item{cop}{outer_nacopula to be estimated}
  \item{interval}{bivariate vector denoting the interval where
    optimization takes place}
  \item{method}{a character string specifying the method (currently only
    for \code{etau()}).}
  \item{\dots}{possibly additional arguments for
    \code{\link{cor}} (\code{etau}), for \code{\link{optimize}}
    (\code{edmle}), or \code{\link{safeUroot}} (\code{ebeta}).}
}
\value{
  \describe{
    \item{ebeta}{Blomqvist beta estimator, the return value of
      \code{\link{safeUroot}}, i.e., almost the same as the value of
      \code{\link{uniroot}}.}
    \item{edmle}{diagonal maximum likelihood estimator, the return value of
      \code{\link{optimize}()}.}
    \item{etau}{averaged pairwise cor() estimators.}
  }
}
\details{
  For \code{etau()}, the \code{method} argument \code{"tau.mean"}
  indicates that the average of the sample versions of Kendall's tau are
  computed first and then theta is determined;
  \code{"theta.mean"} stands for first computing all pairweise Kendall's tau
  estimators and then returning the mean of these estimators.

  For \code{ebeta()}, a method of moment principle equates the
  Blomqvist's beta with the population version of the Archimedean copula.
%% FIXME?  beta.() is currently *not* exported nor documented, hence we
%% can't say more above.
}
\author{Marius Hofert}
\seealso{
  \code{\link{enacopula}}, also for references and examples;
  Specifically, the more sophisticated estimators,
  \code{\link{emle}} (Maximum Likelihood) and
  \code{\link{emde}} (Minimum Distance).
}
\examples{
%% FIXME
}
\keyword{models}
