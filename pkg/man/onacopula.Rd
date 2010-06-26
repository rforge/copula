\name{onacopula}
\alias{onacopula}
\title{Constructing (Outer) Nested Archimedean Copulas}
\description{
  Constructing (outer) nested Archimedean copulas (class
  \code{\linkS4class{outer_nacopula}}) is most conveniently done via
  \code{onacopula()}, using a nested \eqn{C(...)} notation.
}
\usage{
onacopula(family, nacStructure)
}
\arguments{
  \item{family}{either a \code{\link{character}} string, the short or longer form
    of the Archimedean family name (e.g., "Clayton" or simply "C"), see
    the \code{\link{acopula-families}} documentation,
    \cr
    or an \code{\linkS4class{acopula}} family object.}
  \item{nacStructure}{a \dQuote{formula} of the form
    \deqn{C(\vartheta, c(i_1,\dots,i_c), list(C(..), ..., C(..))).}{C(th, c(i1,..ic), list(C(..), ..., C(..))).}
    Note that \eqn{C()} has (maximally) three arguments: the first is
    the copula parameter (vector) \eqn{\vartheta}{theta}, the second a
    (possibly empty) vector of integer indices of components (for the
    \code{comp} slot in \code{\linkS4class{nacopula}}s), and finally a
    (possibly empty) list of child copulas, each specified with in the
    \eqn{C(..)} notation themselves.}
}
\value{
  An outer nested Archimedean copula object, i.e., of class
  \code{\linkS4class{"outer_nacopula"}}.
}
\author{Martin Maechler}
\seealso{
  The class definitions \code{\linkS4class{"nacopula"}},
  \code{\linkS4class{"outer_nacopula"}}, and \code{\linkS4class{"acopula"}}.
}
\references{
	Hofert, M. (2010a),
	Efficiently sampling nested Archimedean copulas,
	\emph{Computational Statistics & Data Analysis}, in press.

	Hofert, M. (2010b),
	\emph{Sampling Nested Archimedean Copulas with Applications to CDO Pricing},
	Suedwestdeutscher Verlag fuer Hochschulschriften AG & Co. KG.
}
\examples{
## Construct a three-dimensional nested Gumbel copula with parameters chosen
## such that the Kendall's tau of the respective bivariate margins are
## 0.2 and 0.5.
theta0 <- copGumbel@tauInv(.2)
theta1 <- copGumbel@tauInv(.5)
C3 <- onacopula("G", C(theta0, 1, C(theta1, c(2,3))))

## Good error message if the component ("coordinate") indices are wrong / do not match:
err <- try(onacopula("G", C(theta0, 2, C(theta1, c(3,2)))))

## Compute the probability of falling in the cube [0,.01]^3 for this copula
pnacopula(C3, rep(.01,3))

## Compute the probability of falling in the cube [.99,1]^3
prob(C3, rep(.99, 3), rep(1, 3))
}
\keyword{multivariate}
