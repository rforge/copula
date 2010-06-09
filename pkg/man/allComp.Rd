\name{allComp}
\alias{allComp}
\title{All Components of a (Inner or Outer) Nested Archimedean Copula}
\description{
  Return Integer Vector of the \emph{indices} of all components for
  which the nested archimedean copula \code{x} is defined.
  %% TODO, the above is too imprecise
  \cr
  This is typically only used by programmers investigating the exact
  nesting structure.

  For an \code{\linkS4class{outer_nacopula}} object
  \code{x}, \code{allComp(x)} must be the same as
  \code{1:\link{dim}(x)}, whereas its \dQuote{inner} component copulas
  will each contain a \emph{subset} of those indices only.
}
\usage{
allComp(x)
}
\arguments{
  \item{x}{an \R object of inheriting from class \code{\linkS4class{nacopula}}.}
}
\value{
  an integer vector of indices \eqn{j} of all components \eqn{u_j} to
  which the copula applies; .
}
\author{Martin Maechler}
\examples{
 example(onacopula)
 allComp(C3d.AMH) # -> 1:3
 allComp(C3d.AMH@childCops[[1]]) # only  (2, 3)
}
\keyword{manip}
\keyword{utilities}
