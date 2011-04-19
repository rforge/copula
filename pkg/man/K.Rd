\name{K}
\alias{K}
\title{Kendall distribution function}
\description{
  Kendall distribution function
%% FIXME: Say more ...
}
\usage{
K(t, cop, d, n.MC = 0)
}
\arguments{
  \item{t}{evaluation point(s)}
  \item{cop}{acopula with specified parameter}
  \item{d}{dimension}
  \item{n.MC}{if > 0 a Monte Carlo approach is applied with sample size equal
    to n.MC; otherwise the exact formula is used.}
}
\value{Kendall distribution function at t}
\author{Marius Hofert}
% \details{
% ...
% }
\seealso{
  the copula families, e.g., \code{\link{copGumbel}}, also for references.
}
\examples{
## TODO
}
\keyword{distribution}
