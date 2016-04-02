\name{pairs2}
\alias{pairs2}
\title{Scatter-Plot Matrix (pairs) with Nice Defaults}
\description{
  A version of \pkg{graphics}' \code{\link{pairs}()} function,
  particularly useful for visualizing dependence in multivariate
  (copula) data.
}
\usage{
pairs2(x, labels = NULL, labels.null.lab = "U", row1attop = FALSE,
       pch = ".", \dots)
}
\arguments{
  \item{x}{A numeric matrix or as.matrix(.)able.}
  \item{labels}{The variable names, typically unspecified.}
  \item{labels.null.lab}{The \code{\link{character}} string
    determining the \dQuote{base name} of the variable labels in case
    \code{labels} is \code{NULL} and \code{x} does not have all column names given.}
  \item{row1attop}{See \code{\link{pairs}()}.}
  \item{pch}{The plot character used.}
  \item{\dots}{Further arguments, passed to \code{\link{pairs}()}.}
}
\value{\code{\link{invisible}()}
}
%\author{Marius Hofert}
\seealso{
  \code{\link{splom2}()} for a similar function based on
  \code{\link[lattice]{splom}()}.
}
\examples{
## Create a 100 x 7 matrix of random variates from a t distribution
## with four degrees of freedom and plot the generated data
U <- matrix(rt(700, df = 4), ncol = 7)
pairs2(U)
}
\keyword{hplot}
