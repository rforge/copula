\name{splom2}
\alias{splom2}
\title{Scatter-Plot Matrix (splom) with Nice Defaults}
\description{
  A version of \pkg{lattice}'s \code{\link{splom}()} function,
  particularly useful for visualizing dependence in multivariate
  (copula) data.
}
\usage{
splom2(x, varnames = NULL, varnames.null.lab = "U",
       xlab = "", col.mat = NULL, bg.col.mat = NULL, \dots)
}
\arguments{
  \item{x}{A numeric matrix or as.matrix(.)able.}
  \item{varnames}{The variable names, typically unspecified.}
  \item{varnames.null.lab}{The \code{\link{character}} string
    determining the \dQuote{base name} of the variable labels in case
    \code{varnames} is \code{NULL} and \code{x} does not have all column names given.}
  \item{xlab}{The x-axis label.}
  \item{col.mat}{A matrix of colors for the plot symbols (the default is the
    setting as obtained from
    \code{\link[lattice]{trellis.par.get}("plot.symbol")$col}).}
  \item{bg.col.mat}{A matrix of colors for the background (the default is the
    setting as obtained from
    \code{\link[lattice]{trellis.par.get}("background")$col}).}
  \item{\dots}{Further arguments, passed to \code{\link[lattice]{splom}()}.}
}
\value{From \code{\link{splom}()}, an \R object of class \code{"trellis"}.
}
\author{Martin Maechler, with a hint from Deepayan Sarkar; based on
  ideas from Marius Hofert.}
\seealso{
  \code{\link{pairs2}()} for a similar function based on
  \code{\link[graphics]{pairs}()}.
}
\examples{
## Create a 100 x 7 matrix of random variates from a t distribution
## with four degrees of freedom and plot the generated data
U <- matrix(rt(700, df = 4), ncol = 7)
G <- splom2(U)
G
}
\keyword{hplot}
