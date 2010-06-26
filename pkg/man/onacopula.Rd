\name{onacopula}
\alias{onacopula}
\alias{onacopulaL}
\title{Constructing (Outer) Nested Archimedean Copulas}
\description{
  Constructing (outer) nested Archimedean copulas (class
  \code{\linkS4class{outer_nacopula}}) is most conveniently done via
  \code{onacopula()}, using a nested \eqn{C(...)} notation.

  Slightly less conveniently, but with the option to pass a \code{\link{list}}
  structure, \code{onacopulaL()} can be used, typically from inside
  another function programmatically.
}
\usage{
onacopula (family, nacStructure)
onacopulaL(family, nacList)
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
  \item{nacList}{a \code{\link{list}} of length 3 (or 2), with elements
    \enumerate{
      \item theta: \eqn{\vartheta}
      \item \code{comps}: components \eqn{c(i_1,\dots,i_c)}{c(i1,..ic)}
      \item \code{children}: a \link{list} which must be a
      \code{nacList} itself and may be missing to denote the empty
      \code{list()}.
    }
  }
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

## Equivalent construction with onacopulaL():
str(NAlis <- list(theta0, 1, list(list(theta1, c(2,3)))))
C3. <- onacopulaL("Gumbel", NAlis)
stopifnot(identical(C3, C3.))
}
\keyword{multivariate}
