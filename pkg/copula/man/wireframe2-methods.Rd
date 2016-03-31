\name{wireframe2-methods}
\alias{wireframe2}
\alias{wireframe2-methods}
\alias{wireframe2,matrix-method}
\alias{wireframe2,data.frame-method}
\alias{wireframe2,copula-method}
\alias{wireframe2,mvdc-method}
\docType{methods}
\title{Methods `wireframe2' in Package `copula'}
\description{
  Methods \code{wireframe2()} to draw wireframe plots (of
  two dimensional distributions from package \pkg{copula}).
}
\usage{
\S4method{wireframe2}{matrix}(x, labels = NULL,
      labels.null.lab = parse(text = c("u[1]", "u[2]", "C(u[1],u[2])")),
      alpha.regions = 0.5, scales = list(arrows = FALSE, col = 1),
      par.settings = standard.theme(color = FALSE), \dots)
\S4method{wireframe2}{data.frame}(x, labels = NULL,
      labels.null.lab = parse(text = c("u[1]", "u[2]", "C(u[1],u[2])")),
      alpha.regions = 0.5, scales = list(arrows = FALSE, col = 1),
      par.settings = standard.theme(color = FALSE), \dots)
\S4method{wireframe2}{copula}(x, FUN, n.grid = 26, delta = 0,
      labels = c("u[1]", "u[2]", deparse(substitute(FUN))[1]), \dots)
\S4method{wireframe2}{mvdc}(x, FUN, xlim, ylim, n.grid = 26,
      labels = c("x[1]", "x[2]", deparse(substitute(FUN))[1]), \dots)
}
\arguments{
  \item{x}{A \code{"\linkS4class{matrix}"},
    \code{"\linkS4class{data.frame}"}, \code{"\linkS4class{copula}"} or
    a \code{"\linkS4class{mvdc}"} object.}
  \item{labels}{The x-, y- and z-axis labels.}
  \item{labels.null.lab}{A character vector of length 3 determining the
    \dQuote{base name} of the variable labels in case \code{x} does not
    have all column names given and \code{labels} is \code{NULL}.}
  \item{alpha.regions}{See \code{\link{wireframe}()}.}
  \item{scales}{See \code{\link{wireframe}()}.}
  \item{par.settings}{See \code{\link{wireframe}()}.}
  \item{FUN}{The \code{\link{function}} to be plotted; typically
    \code{\link{dCopula}} or \code{\link{pCopula}}.}
  \item{n.grid}{The number of grid points used in each dimension. This
    can be a vector of length two, giving the number of grid points used
    in x- and y-direction, respectively; the function \code{FUN} will be
    evaluated on the corresponding (x,y)-grid.}
  \item{delta}{A small number in \eqn{[0, \frac 1 2)}{[0, 1/2)}
    influencing the evaluation boundaries. The x- and y- vectors will
    have the range \code{[0+delta, 1-delta]}, the default being \code{[0,1]}.}
  \item{xlim, ylim}{The \code{\link{range}} of the x and
    y variables, respectively.}
  \item{\dots}{Additional arguments passed to the underlying
    \code{\link{wireframe}()}.}
}
\section{Methods}{
  Wireframe plots for objects of class \code{"\linkS4class{matrix}"} ,
  \code{"\linkS4class{data.frame}"}, \code{"\linkS4class{copula}"} or
  \code{"\linkS4class{mvdc}"}.
}
\value{
  An object of class \dQuote{trellis} as returned by \code{\link[lattice]{wireframe}()}.
}
\author{Marius Hofert}
\seealso{
  The \code{\link{persp-methods}} for drawing perspective plots via base graphics.
}
\examples{
## For 'matrix' objects
## 1) The independence copula Pi
n.grid <- 26
u <- seq(0, 1, length.out = n.grid)
grid <- expand.grid("u[1]" = u, "u[2]" = u)
x <- cbind(grid, "Pi(u[1],u[2])" = apply(grid, 1, prod)) # evaluate Pi on 'grid'
wireframe2(x) # labels are parsed

## 2) The Frechet--Hoeffding bounds W and M
W <- function(u) pmax(0, rowSums(u)-1) # lower bound W
M <- function(u) apply(u, 1, min) # upper bound M
x.W <- cbind(grid, "W(u[1],u[2])" = W(grid)) # evaluate W on 'grid'
x.M <- cbind(grid, "M(u[1],u[2])" = M(grid)) # evaluate M on 'grid'
wireframe2(x.W) # plot of W
wireframe2(x.M) # plot of M

## For 'copula' objects
cop <- frankCopula(-4)
wireframe2(cop, pCopula)
wireframe2(cop, dCopula)
wireframe2(cop, dCopula, delta=0.01)
wireframe2(cop, pCopula, labels=c("x", "y", "z"))

## For 'mvdc' objects
mvNN <- mvdc(gumbelCopula(3), c("norm", "norm"),
             list(list(mean = 0, sd = 1), list(mean = 1)))
wireframe2(mvNN, dMvdc, xlim=c(-2, 2), ylim=c(-1, 3))
}
\keyword{methods}
\keyword{hplot}
