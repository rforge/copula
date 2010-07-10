\name{splom2}
\alias{splom2}
\title{Scatterplot Matrix [SPLOM] With Nice Variable Names}
\description{
  A version of \pkg{lattice}'s \code{\link{splom}} function,
  particularly useful for visualizing multivariate data sampled from
  copulas, notably nested Archimedean ones.

  \bold{Experimental} currently!  We may replace the interface entirely,
  e.g., to accept an \code{"\linkS4class{outer_nacopula}"} object from
  which to sample.
}
\usage{
splom2(data, varnames, Vname="U", ...)
}
\arguments{
  \item{data}{numeric matrix or as.matrix(.)able}
  \item{varnames}{variable names, typically unspecified}
  \item{Vname}{character string to become the "base name" of the variables}
  \item{...}{further arguments to splom()}
}
\value{an splom() object}
\author{Martin Maechler, with a hint from Deepayan Sarkar; based on
  ideas from Marius Hofert.}
\examples{
## Create a 100 x 7 matrix of random variates from a t distribution
## with four degrees of freedom and plot the generated data
U7 <- matrix(rt(700, 4), 100, 7)
G <- splom2(U7)
G
}
\keyword{hplot}
