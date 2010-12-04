\name{polynEval}
\alias{polynEval}
\title{Evaluate Polynomials}
\description{
  Evaluate a univariate polynomial at \code{x},
  i.e., compute \code{coef[1] + coef[2]*x + ... + coef[p+1]* x^p},
  where \code{x} is typically a vector.
}
\usage{
polynEval(coef, x)
}
\arguments{
  \item{coef}{numeric vector.  If a vector, \code{x} can be an
    array and the result matches \code{x}.}
  \item{x}{numeric vector or array.}
}
\details{
  The stable \dQuote{Horner rule} is used for evaluation.

  Using the C code speeds up the already fast \R code available in
  \code{\link[sfsmisc]{polyn.eval}()} in package \pkg{sfsmisc}.
}
\value{
  numeric vector or array, with the same dimensions as \code{x},
  containing the polynomial values \eqn{p(x)}.
}
\author{Martin Maechler; the \R version has been in package
  \pkg{sfsmisc} for ages.}
\seealso{For much more sophisticated handling of polynomials, use the
  \code{polynom} package, e.g. \code{\link[polynom]{predict.polynomial}}.
}
\examples{
polynEval(c(1,-2,1), x = -2:7)# (x - 1)^2
polynEval(c(0, 24, -50, 35, -10, 1),
          x = matrix(0:5, 2,3))# 5 zeros!
}
\keyword{arith}
