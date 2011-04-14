\name{emde}
\alias{emde}
\title{Minimum Distance Estimation for Nested Archimedean Copulas}
\description{
  Minimum distance estimation for nested Archimedean copulas.
  %% FIXME
}
\usage{
emde(u, cop,
     method=c("mde.normal.CvM", "mde.normal.KS",
              "mde.log.CvM",    "mde.log.KS"),
     interval=paraOptInterval(u, cop@copula@name),
     include.K = ncol(u) <= 5, \dots)
}
\arguments{
  \item{u}{numeric matrix of realizations following the copula}
  \item{cop}{outer_nacopula to be estimated}
  \item{interval}{bivariate vector denoting the interval where
    optimization takes place}
  \item{include.K}{logical indicating whether the last component, K, is
    also used or not.}
  \item{method}{a \code{\link{character}} string specifying the distance
    method. One (or a unique abbreviation) of one of
    \describe{
      \item{\code{"mde.normal.CvM"}}{map to a chi-square distribution
	(Cramer-von Mises distance).}
      \item{\code{"mde.normal.KS"}}{map to a chi-square distribution
	(Kolmogorov-Smirnov distance).}
      \item{\code{"mde.log.CvM"}}{map to an Erlang distribution
	(Cramer-von Mises distance).}
      \item{\code{"mde.log.KS"}}{map to an Erlang distribution
	(Kolmogorov-Smirnov distance).}
    }
  }
  \item{\dots}{possibly additional arguments for
    \code{\link{optimize}}.}
}
\value{
  the minimum distance estimator, specifically the \code{\link{list}} as
  returned by \code{\link{optimize}()}.
}
\author{Marius Hofert}
\seealso{
  \code{\link{enacopula}}, also for references and examples.
}
\keyword{models}
