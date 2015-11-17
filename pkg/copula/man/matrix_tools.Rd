\name{matrix_tools}
\alias{p2P}
\alias{P2p}
\alias{getSigma}
\alias{extreme_pairs}
\title{Tools to Work with Matrices}
\description{
  \code{p2P()} creates a \code{\link{matrix}} from a given
  \code{\link{vector}} of parameters.  \code{P2p()} creates a
  numeric vector from a given \code{\link{matrix}},
  currently useful for elliptical copulas.

  \code{getSigma()} returns the \eqn{d \times d}{d * d} symmetric
  matrix \eqn{\Sigma} which is called \dQuote{Rho} as well, written
  (capital Greek \eqn{\rho}{rho} !) as \eqn{P} (and hence sometimes
  erronously pronounced "Pee").
  Note that \code{getSigma()} works for all elliptical copulas and uses
  \code{p2P()} for the \dQuote{unstuctured} case, \code{dispstr = "un"}.

  \code{extreme_pairs()} identifies pairs with the smallest (or largest,
  if \code{decreasing=TRUE}) \code{n} values in a symmetric matrix.
}
\usage{
p2P(param, d)
P2p(P)
getSigma(copula)
extreme_pairs(x, n=6, decreasing=FALSE)
}
\arguments{
  \item{param}{A parameter vector.}
  \item{d}{Dimension of the resulting \code{\link{matrix}}.}
  \item{P}{A \code{\link{matrix}} which should be converted to a vector.}
  \item{copula}{An \bold{elliptical} copula, i.e., an object (extending)
    class \code{\linkS4class{ellipCopula}}; typically resulting from
    \code{\link{tCopula}()} or \code{\link{normalCopula}()}.}
  \item{x}{A symmetric \code{\link{matrix}}.}
  \item{n}{The number of pairs with smallest (or largest) values to be
    displayed.}
  \item{decreasing}{A \code{\link{logical}} indicating whether pairs
    corresponding to the smallest (for \code{decreasing=FALSE}; the
    default) or largest (for \code{decreasing=TRUE}) \code{n} values
    should be displayed.}
}
\value{
    \describe{
    \item{\code{p2P}:}{A symmetric \code{\link{matrix}} with ones on the
      diagonal and the values of \code{param} filled column-wise below the
      diagonal (which corresponds to row-wise filling above the diagonal).}
    \item{\code{P2p}:}{\code{\link{vector}} of column-wise below-diagonal
      entries of \code{P} (equal to the row-wise above-diagonal entries in
      case of a symmetric matrix).}
    \item{\code{getSigma}:}{\code{\link{matrix}} as from \code{p2P()}
      for all cases of elliptical copulas.}
    \item{\code{extreme_pairs}:}{A (\code{n}, 3)-\code{\link{matrix}}
      containing the row and column indices (first two columns) of the pairs
      with the most extreme (smallest or largest) entries in \code{x}.}
  }
}
\details{
  These auxiliary functions are often used when working with
  elliptical copulas.
}
\seealso{
  \code{\link{ellipCopula}}, \code{\link{tCopula}}, \code{\link{normalCopula}}.
}
%\author{Marius Hofert.}
\examples{
## display the definitions
p2P
P2p
extreme_pairs

param <- (2:7)/10
tC <- tCopula(param, dim = 4, dispstr="un", df = 3)
## consistency of the three functions :
P <- p2P(param, d=4)
stopifnot(identical(param, P2p(P)),
	  identical(P, getSigma(tC)))

## Toeplitz case:
(tCt <- tCopula((2:6)/10, dim = 6, disp = "toep"))
(rhoP <- tCt@getRho(tCt))
stopifnot(identical(getSigma (tCt),
		    toeplitz (c(1, rhoP))))
}
\keyword{array}
\keyword{manip}

