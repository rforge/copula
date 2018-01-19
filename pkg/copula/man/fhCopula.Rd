\name{fhCopula}
\title{Construction of Fr\enc{é}{e}chet-Hoeffding Bound Copula Objects}
\alias{fhCopula}
\alias{lowfhCopula}
\alias{upfhCopula}
\description{
  Constructs the Fr\enc{é}{e}chet-Hoeffding bound copulas.
}
\usage{
fhCopula(family = c("upper", "lower"), dim = 2L)

lowfhCopula(dim = 2L)
 upfhCopula(dim = 2L)
}
\arguments{
  \item{family}{a character string specifying the
    Fr\enc{é}{e}chet-Hoeffding bound copula.}
  \item{dim}{the dimension of the copula; note that the lower
    Fr\enc{é}{e}chet-Hoeffding bound is only available in the bivariate case.}
}
\value{
  A copula object of class \code{"\linkS4class{lowfhCopula}"}
  or \code{"\linkS4class{upfhCopula}"}.
}
\note{
  \code{fhCopula()} is a wrapper for \code{lowfhCopula()} and
  \code{upfhCopula()}.

  The \code{\link{dCopula}()} method will simply return an error
  for these copulas (as their density does not exist). Also,
  since the Fr\enc{é}{e}chet-Hoeffding bound copulas are not
  parametric, certain methods available for parametric copulas are
  not available.
}
%\author{Marius Hofert}
\examples{
try(W <- lowfhCopula(dim = 3))
W <- lowfhCopula()
wireframe2(W, FUN = pCopula)

M <- upfhCopula(dim = 2)
wireframe2(M, FUN = pCopula)

M <- upfhCopula(dim = 3)
set.seed(271)
splom2(M, n = 100)
}
\keyword{distribution}
\keyword{multivariate}
