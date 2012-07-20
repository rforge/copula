\name{archmCopula}
\alias{archmCopula}
\alias{claytonCopula}
\alias{frankCopula}
\alias{gumbelCopula}
\alias{amhCopula}
\title{Construction of Archimedean copula class object}
\description{
  Constructs an Archimedean copula class object with its corresponding
  parameter and dimension.
}
\usage{
archmCopula(family, param = NA_real_, dim = 2, ...)

claytonCopula(param = NA_real_, dim = 2)
frankCopula(param = NA_real_, dim = 2)
gumbelCopula(param = NA_real_, dim = 2)
amhCopula(param = NA_real_, dim = 2)
}
\arguments{
  \item{family}{a character string specifying the family of an
    Archimedean copula.  Currently supported families are "clayton",
    "frank", "amh", and "gumbel".
  }
  \item{param}{number (\code{\link{numeric}}) specifying the copula parameter.}
  \item{dim}{the dimension of the copula.}
  \item{\dots}{ currently nothing. }
}
\value{
  An Archimedean copula object of class \code{"claytonCopula"},
  \code{"frankCopula"}, \code{"gumbelCopula"}, or \code{"amhCopula"}.
}
\details{
  \code{archmCopula()} is a wrapper for \code{claytonCopula()},
  \code{frankCopula()}, \code{gumbelCopula()} and \code{amhCopula()}.

  Only the bivariate Ali-Mikhail-Haq copula family (\code{"amhCopula"})
  is available.

  The maximum dimension for which the expression of the pdf is available
  is 6 for the Clayton, Gumbel and Frank families. The cdf expression is
  always available.

  The maximum dimension for which \code{"dcopula"} can be evaluated is
  10 for the Clayton and Gumbel families, and 6 for the Frank family.
  They are also the maximum dimensions for which maximum likelihood
  estimation can be done.

  It is now recommended to work with the
  \code{\linkS4class{acopula}}-classed Archimedean copulas, as there is
  now restriction on the dimension there.
}
\references{
  R.B. Nelsen (2006), \emph{An introduction to Copulas}, Springer, New York.
}
\seealso{
  \code{\link{ellipCopula}}, \code{\link{evCopula}}.

  \code{\linkS4class{acopula}}-classed Archimedean copulas, such as
  \code{\link{copClayton}}, \code{\link{copGumbel}}, etc.
}
\examples{
clayton.cop <- claytonCopula(2, dim = 3)
## scatterplot3d(rcopula(clayton.cop, 1000))

frank.cop <- frankCopula(3)
persp(frank.cop, dcopula)

gumbel.cop <- archmCopula("gumbel", 5)
contour(gumbel.cop, dcopula)

amh.cop <- amhCopula(0.5)
u. <- as.matrix(expand.grid(u=(0:10)/10, v=(0:10)/10, KEEP.OUT.ATTRS=FALSE))
du <- dcopula(amh.cop, u.)
stopifnot(is.finite(du) | apply(u. == 0, 1,any)| apply(u. == 1, 1,any))

## A 7-dim Frank copula
frank.cop <- frankCopula(3, dim = 7)
x <- rcopula(frank.cop, 5)
## dcopula now *does* work:
dcopula(frank.cop, x)

## A 7-dim Gumbel copula
gumbel.cop <- gumbelCopula(2, dim = 7)
dcopula(gumbel.cop, x)
}
\keyword{distribution}
\keyword{multivariate}
