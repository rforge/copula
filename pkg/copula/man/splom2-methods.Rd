\name{splom2-methods}
\alias{splom2}
\alias{splom2-methods}
\alias{splom2,matrix-method}
\alias{splom2,data.frame-method}
\alias{splom2,copula-method}
\alias{splom2,mvdc-method}
\docType{methods}
\title{Methods `splom2' in Package `copula'}
\description{
  Methods \code{splom2()} to draw scatter-plot matrices of
  distributions from package \pkg{copula}.
}
\usage{
\S4method{splom2}{matrix}(x, varnames = NULL, varnames.null.lab = "U",
      xlab = "", col.mat = NULL, bg.col.mat = NULL, \dots)
\S4method{splom2}{data.frame}(x, varnames = NULL, varnames.null.lab = "U",
      xlab = "", col.mat = NULL, bg.col.mat = NULL, \dots)
\S4method{splom2}{copula}(x, n, \dots)
\S4method{splom2}{mvdc}(x, n, varnames.null.lab = "X", \dots)
}
\arguments{
  \item{x}{A \code{"\linkS4class{matrix}"},
    \code{"\linkS4class{data.frame}"}, \code{"\linkS4class{copula}"} or
    a \code{"\linkS4class{mvdc}"} object.}
  \item{n}{The sample size.}
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
  \item{\dots}{Additional arguments passed to the underlying
    \code{\link[lattice]{splom}()}.}
}
\value{From \code{\link{splom}()}, an \R object of class \code{"trellis"}.
}
%\author{Marius Hofert and Martin Maechler, with a hints from Deepayan Sarkar.}
\seealso{
  \code{\link{pairs2}()} for a similar function (for matrices and data
  frames) based on \code{\link[graphics]{pairs}()}.
}
\examples{
## For 'matrix' objects
## Create a 100 x 7 matrix of random variates from a t distribution
## with four degrees of freedom and plot the generated data
n <- 1000 # sample size
d <- 3 # dimension
nu <- 4 # degrees of freedom
tau <- 0.5 # Kendall's tau
th <- iTau(tCopula(df = nu), tau) # corresponding parameter
cop <- tCopula(th, dim = d, df = nu) # define copula object
set.seed(271)
U <- rCopula(n, copula = cop)
splom2(U)

## For 'copula' objects
set.seed(271)
splom2(cop, n = n) # same as above

## For 'mvdc' objects
mvNN <- mvdc(cop, c("norm", "norm", "exp"),
             list(list(mean = 0, sd = 1), list(mean = 1), list(rate = 2)))
splom2(mvNN, n = n)
}
\keyword{hplot}
