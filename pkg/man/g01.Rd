\name{g01}
\alias{g01}
\alias{g01trafo}
\title{Transformation and Anderson-Darling Test of Multivariate Uniform U[0,1]}
\description{
  \code{g01()} performs the Anderson-Darling test for U[0,1] of
  supposedly U[0,1]^d distributed data after that is transformed to
  univariate via \code{g01trafo()}.

  \code{g01trafo()} transforms supposedly \eqn{U[0,1]^d} distributed vectors of
  random variates to \eqn{U[0,1]}--distributed variates (for checking
  uniformity in a one-dimensional setup).
}
\usage{
g01     (u, method = c("normal", "log"))
g01trafo(u, method = c("normal", "log"))
}
\arguments{
  \item{u}{matrix, say \eqn{n\times d}{n * d}, of random variates in
    \eqn{[0,1]} to be transformed or tested, respectively.}
  \item{method}{character string; either
    \code{"normal"}, to map to a chi-square distribution via cumulative
    normal transform (\code{\link{pnorm}}) or \cr
    \code{"log"}, to map to an Erlang distribution via \code{-log(u)}.
  }
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  \code{g01()} returns the Anderson-Darling test result, i.e., the value
  of \code{\link[ADGofTest]{ad.test}()} from package \pkg{ADGofTest}.

  \code{g01trafo()} returns a numeric vector of length \eqn{n}
  (=\code{nrow(u)}) of 1D transformed values.
}
\references{
  ...%% FIXME
}
\seealso{
  ...%% FIXME
}
\examples{
## FIXME -- add some
}
\keyword{distribution}
\keyword{multivariate}
