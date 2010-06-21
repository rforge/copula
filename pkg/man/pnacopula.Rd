\name{pnacopula}
\docType{methods}
\alias{pnacopula}
\title{Evaluation of Nested Archimedean Copulas}
\description{
  \code{pnacopula} evaluates a nested Archimedean copula (object of class
  \code{\linkS4class{nacopula}}) at the given vector \code{u}.
}
\usage{
pnacopula(x, u)
}
\arguments{
  \item{x}{a (nested) Archimedean copula of dimension \eqn{d}, i.e., an
    object of class \code{\linkS4class{nacopula}}, typically from
    \code{\link{onacopula}(..)}.}
  \item{u}{a \code{\link{numeric}} vector of length \eqn{d}.}
}
\value{
  A \code{\link{numeric}} in \eqn{[0,1]} which is the copula evaluated
  at \code{u}. (Currently not parallelized.)
}
\author{Marius Hofert, Martin Maechler}
\examples{
## Construct a three-dimensional nested Joe copula with parameters
## chosen such that the Kendall's tau of the respective bivariate margins
## are 0.2 and 0.5.
theta0 <- copJoe@tauInv(.2)
theta1 <- copJoe@tauInv(.5)
C3 <- onacopula("J", C(theta0, 1, C(theta1, c(2,3))))

## Evaluate this copula at the point u
u <- c(.7,.8,.6)
pnacopula(C3, u)

## Evaluate this copula at several points
u <- matrix(runif(300), ncol=3)
apply(u, 1, pnacopula, x=C3)
}
\keyword{methods}
