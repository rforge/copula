\name{p2P}
\alias{p2P}
\alias{P2p}
\title{Auxiliary Functions to Convert Matrices of Parameters to Vectors
  and Vice Versa}
\description{
  \code{p2P()} creates a \code{\link{matrix}} from a given
  \code{\link{vector}} of parameters. \code{P2p()} creates a
  \code{\link{vector}} from a given \code{\link{matrix}}.
}
\usage{
p2P(param, d)
P2p(P)
}
\arguments{
  \item{param}{a parameter vector.}
  \item{d}{dimension of the resulting \code{\link{matrix}}.}
  \item{P}{\code{\link{matrix}} which should be converted to a vector.}
}
\value{
  \describe{
    \item{\code{p2P}: }{a symmetric \code{\link{matrix}} with ones on the
      diagonal and the values of \code{param} filled column-wise below the
      diagonal (which corresponds to row-wise filling above the diagonal).}
    \item{\code{P2p}: }{\code{\link{vector}} of column-wise below-diagonal
      entries of \code{P} (equal to the row-wise above-diagonal entries in
      case of a symmetric matrix).}
  }
}
\details{
  These two auxiliary functions are often used when working with
  elliptical copulas.
}
\author{Marius Hofert.}
\examples{
## display the definitions
p2P
P2p
## nice check
param <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
stopifnot(P2p(p2P(param, 4))==param)
}
\keyword{parameter}
\keyword{vector}
\keyword{matrix}
