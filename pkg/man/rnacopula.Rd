\name{rnacopula}
\alias{rnacopula}
\title{Sampling Nested Archimedean Copulas}
\description{
  Random number generation from nested Archimedean copulas (of class
  \code{\linkS4class{outer_nacopula}}, specifically), aka
  \dQuote{sampling} nested Archimedean copulas will generate \code{n}
  random vectors of dimension \eqn{d} (= \code{dim(x)}).
}
\usage{
rnacopula(n, x, ...)
}
\arguments{
  \item{x}{an \R object of \code{\link{class}}
    \code{"\linkS4class{outer_nacopula}"}, typically from \code{\link{onacopula}()}.}
  \item{n}{integer specifying the sample size, i.e., the number of
    copula-distributed random vectors \eqn{\mathbf{U}_i}{U_i}, to be generated.}
  \item{\dots}{possibly further arguments for the given copula family.}
}
\details{
  The generation happens by calling \code{\link{rnchild}()} on
  each child copula (which itself recursively descends the tree implied
  by the nested Archimedean structure). The algorithm is based on a
  mixture representation of the generic distribution functions
  \eqn{F_{0}}{F0} and \eqn{F_{01}}{F01} and is presented in
  McNeil~(2008). Details about how to efficiently sample the distribution
  functions \eqn{F_{0}}{F0} and \eqn{F_{01}}{F01} can be found in
  Hofert~(2010a) and Hofert (2010b).
}
\value{
  a \code{\link{numeric}} matrix containing the generated vectors of
  random variates from the nested Archimedean copula object \code{x}.
}
\author{Marius Hofert, Martin Maechler}
\seealso{
  \code{\link{rnchild}}; classes \code{"\linkS4class{nacopula}"} and
  \code{"\linkS4class{outer_nacopula}"}, see also \code{\link{onacopula}()}.
}
\references{
  McNeil, A. J. (2008),
  Sampling nested Archimedean copulas,
  \emph{Journal of Statistical Computation and Simulation} \bold{78}, 6, 567--581.

  Further, those of the Archimedean families, e.g., \code{\link{copGumbel}}.
}
\examples{
## Construct a three-dimensional nested Clayton copula with parameters
## chosen such that the Kendall's tau of the respective bivariate margins
## are 0.2 and 0.5 :
C3 <- onacopula("C", C(copClayton@tauInv(0.2), 1,
                       C(copClayton@tauInv(0.5), c(2,3))))
C3

## Sample n vectors of random variates from this copula. This involves
## sampling exponentially tilted stable distributions
n <- 1000
U <- rnacopula(n, C3)

## Plot the drawn vectors of random variates
splom2(U)
}
\keyword{distribution}

