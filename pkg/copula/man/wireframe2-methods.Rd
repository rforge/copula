\name{wireframe2-methods}
\alias{wireframe2-methods}
\alias{wireframe2,copula-method}
\alias{wireframe2,mvdc-method}
\docType{methods}
\title{Methods for Function `wireframe2' in Package `copula'}
\description{
  Methods for function \code{\link{wireframe2}} to draw wireframe plots (of
  two dimensional distributions from package \pkg{copula}).
}
\usage{
\S4method{wireframe2}{copula}(x, FUN, n.grid = 26, delta = 0,
      labels = c("u[1]", "u[2]", deparse(substitute(FUN))[1]), \dots)
\S4method{wireframe2}{mvdc}(x, FUN, xlim, ylim, n.grid = 26,
      labels = c("x[1]", "x[2]", deparse(substitute(FUN))[1]), \dots)
}
\arguments{
  \item{x}{A \code{"\linkS4class{copula}"} or
    a \code{"\linkS4class{mvdc}"} object.}
  \item{FUN}{The \code{\link{function}} to be plotted; typically
    \code{\link{dCopula}} or \code{\link{pCopula}}.}
  \item{n.grid}{The number of grid points used in each dimension. This
    can be a vector of length two, giving the number of grid points used
    in x- and y-direction, respectively; the function \code{FUN} will be
    evaluated on the corresponding (x,y)-grid.}
  \item{delta}{A small number in \eqn{[0, \frac 1 2)}{[0, 1/2)}
    influencing the evaluation boundaries. The x- and y- vectors will
    have the range \code{[0+delta, 1-delta]}, the default being \code{[0,1]}.}
  \item{labels}{A vector of length 3 giving the axis labels.}
  \item{xlim, ylim}{The \code{\link{range}} of the x and
    y variables, respectively.}
  \item{\dots}{Additional arguments passed to the underlying
    \code{\link{wireframe2}()}.}
}
\section{Methods}{
  Wireframe plots for both \code{"\linkS4class{copula}"} or
  \code{"\linkS4class{mvdc}"} objects, see \code{x} in the
  \emph{Arguments} section.
}
\value{
  An object of class \dQuote{trellis} as returned by \code{\link[lattice]{wireframe}()}.
}
\author{Marius Hofert}
\seealso{
  The \code{\link{persp-methods}} for drawing perspective plots via base graphics.
}
\examples{
## For 'copula' objects
wireframe2(frankCopula(-0.8), dCopula)
wireframe2(claytonCopula(2),  pCopula, title = "CDF of claytonCopula(2)")

## For 'mvdc' objects
mvNN <- mvdc(gumbelCopula(3), c("norm", "norm"),
          list(list(mean = 0, sd = 1), list(mean = 1)))
wireframe2(mvNN, dMvdc, xlim=c(-2, 2), ylim=c(-1, 3))
}
\keyword{methods}
\keyword{hplot}
