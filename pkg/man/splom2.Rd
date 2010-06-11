\name{splom2}
\alias{splom2}
\title{Scatterplot Matrix [SPLOM] With Nice Variable Names}
\usage{
splom2(data, varnames, Vname="U", ...)
}
\arguments{
  \item{data}{numeric matrix or as.matrix(.)able}
  \item{varnames}{variable names, typically unspecified}
  \item{Vname}{character string to become "root variable name"}
  \item{...}{further arguments to splom()}
}
\value{a splom() object}
\author{Martin Maechler, with a hint from Deepayan Sarkar}
\examples{
U7 <- matrix(rt(700, 4), 100,7)
G <- splom2(U7)
G
}
\keyword{hplot}
