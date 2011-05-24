\name{paraOptInterval}
\alias{paraOptInterval}
\title{Initial Interval for Archimedean Copula Estimation}
\description{
  Compute an initial interval for optimization/estimation routines (only a
  heuristic; if this fails, choose your own interval).
}
\usage{
paraOptInterval(u, family, h=0.15)
}
\arguments{
  \item{u}{matrix of realizations following a copula}
  \item{family}{Archimedean family}
  \item{h}{for enlarging the tau-interval}
}
\details{
  ... FIXME
}
\value{initial interval which can be used for optimization (e.g., for
  \code{\link{emle}}).
}
% \references{
% }
\author{Marius Hofert}
% \seealso{
% }
\examples{
}
