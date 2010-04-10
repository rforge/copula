\name{onACopula}
\alias{onACopula}
\title{(Outer) Nested Archimedean Copula Creation}
\description{
  Creation of (outer) nested archimedean copulas (class
  \code{\linkS4class{outer_nACopula}}) is most conveniently done via
  \code{onACopula()}, using a nested \eqn{C(...)} notation.
}
\usage{
onACopula(family, nACform)
}
\arguments{
  \item{family}{\code{\link{character}} string: the short or longer form
    of copula family name, see the \code{\link{ACopula-families}}
    documentation.}
  \item{nACform}{a \dQuote{formula} of the form
    \deqn{C(\theta, c(i_1,\dots,i_c), list(C(..), ..., C(..))).}{C(th, c(i1,..ic), list(C(..), ..., C(..))).}
    Note that \eqn{C()} has (maximally) three arguments: the first is
    the copula parameter (vector) \eqn{\theta}{theta}, the second a
    vector (possibly empty) of integer component indices (for the
    \code{comp} slot in \code{\linkS4class{nACopula}}s), and the third a
    list (possibly empty) of sub-copulas, each specified with
    \eqn{C(..)} notation themselves.}
}
\value{
  an outer nested archimedean copula object, i.e., of class
  \code{\linkS4class{"outer_nACopula"}}.
}
% \references{
% %% ~put references to the literature/web site here ~
% }
\author{Martin Maechler}
\seealso{
  The class definitions \code{\linkS4class{"nACopula"}}
  (and \code{"outer_nACopula"}) and \code{\linkS4class{"ACopula"}}.
}
\examples{
C3d.AMH <- onACopula("AMH", C(0.7135, 1, C(0.943, 2:3)))
value(C3d.AMH, u = c(1/4, 1/2, 1/4))

## Good error message if the component ("coordinate") indices are wrong:
err <- try(onACopula("AMH", C(0.7135, 2, C(0.943, c(3,2)))))
}
\keyword{multivariate}
