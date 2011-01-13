\name{onacopula}
\alias{onacopula}
\alias{onacopulaL}
\alias{nacopula}
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
%onacopula (family, nacStructure, envir=environment())
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
    \deqn{C(\theta, c(i_1,\dots,i_c), list(C(..), ..., C(..))).}{C(th, c(i1,..ic), list(C(..), ..., C(..))).}
    Note that \eqn{C()} has (maximally) three arguments: the first is
    the copula parameter (vector) \eqn{\theta}{theta}, the second a
    (possibly empty) vector of integer indices of components (for the
    \code{comp} slot in \code{\linkS4class{nacopula}}s), and finally a
    (possibly empty) list of child copulas, each specified with in the
    \eqn{C(..)} notation themselves.}
  % \item{envir}{the \code{\link{environment}} in which the
  %   \code{nacStructure} \dQuote{formula} is evaluated.}
  \item{nacList}{a \code{\link{list}} of length 3 (or 2), with elements
    \enumerate{
      \item theta: \eqn{\theta}
      \item \code{comps}: components \eqn{c(i_1,\dots,i_c)}{c(i1,..ic)}
      \item \code{children}: a \link{list} which must be a
      \code{nacList} itself and may be missing to denote the empty
      \code{list()}.
    }
  }
}
\value{
  An outer nested Archimedean copula object, i.e., of class
  \code{"\linkS4class{outer_nacopula}"}.
}
\author{Martin Maechler}
\seealso{
  The class definitions \code{"\linkS4class{nacopula}"},
  \code{"\linkS4class{outer_nacopula}"}, and \code{"\linkS4class{acopula}"}.
}
\references{
  Hofert, M. (2011)
  Efficiently sampling nested Archimedean copulas,
  \emph{Computational Statistics & Data Analysis} \bold{55}, 57--70.

  Further, those of the Archimedean families, e.g., \code{\link{copGumbel}}.
}
\examples{
## Construct a ten-dimensional Joe copula with parameter such that Kendall's
## tau equals 0.5
theta <- copJoe@tauInv(0.5)
C10 <- onacopula("J",C(theta,1:10))

## Equivalent construction with onacopulaL():
C10. <- onacopulaL("J",list(theta,1:10))
stopifnot(identical(C10, C10.))

## Construct a three-dimensional nested Gumbel copula with parameters such that
## Kendall's tau of the respective bivariate margins are 0.2 and 0.5.
theta0 <- copGumbel@tauInv(.2)
theta1 <- copGumbel@tauInv(.5)
C3 <- onacopula("G", C(theta0, 1, C(theta1, c(2,3))))

## Equivalent construction with onacopulaL():
str(NAlis <- list(theta0, 1, list(list(theta1, c(2,3)))))
C3. <- onacopulaL("Gumbel", NAlis)
stopifnot(identical(C3, C3.))

## Good error message if the component ("coordinate") indices are wrong / do not
## match:
err <- try(onacopula("G", C(theta0, 2, C(theta1, c(3,2)))))

## Compute the probability of falling in the cube [0,.01]^3 for this copula
pnacopula(C3, rep(.01,3))

## Compute the probability of falling in the cube [.99,1]^3
prob(C3, rep(.99, 3), rep(1, 3))
}
\keyword{multivariate}
