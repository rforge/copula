\name{wireframe2-methods}
\title{Perspective Plots - 'wireframe2' in Package 'copula'}
\alias{wireframe2}
\alias{wireframe2-methods}
\alias{wireframe2,matrix-method}
\alias{wireframe2,data.frame-method}
\alias{wireframe2,Copula-method}
\alias{wireframe2,mvdc-method}
\docType{methods}
\description{
  Generic function and methods \code{wireframe2()} to draw
  (\CRANpkg{lattice}) \code{\link[lattice]{wireframe}} (aka
  \dQuote{perspective}) plots of two-dimensional distributions from
  package \pkg{copula}.
}
\usage{
\S4method{wireframe2}{matrix}(x,
      xlim = range(x[,1], finite = TRUE),
      ylim = range(x[,2], finite = TRUE),
      zlim = range(x[,3], finite = TRUE),
      xlab = NULL, ylab = NULL, zlab = NULL,
      alpha.regions = 0.5, scales = list(arrows = FALSE, col = "black"),
      par.settings = standard.theme(color = FALSE), \dots)
\S4method{wireframe2}{data.frame}(x,
      xlim = range(x[,1], finite = TRUE),
      ylim = range(x[,2], finite = TRUE),
      zlim = range(x[,3], finite = TRUE),
      xlab = NULL, ylab = NULL, zlab = NULL,
      alpha.regions = 0.5, scales = list(arrows = FALSE, col = "black"),
      par.settings = standard.theme(color = FALSE), \dots)
\S4method{wireframe2}{Copula}(x, FUN, n.grid = 26, delta = 0,
      xlim = 0:1, ylim = 0:1, zlim = NULL,
      xlab = quote(u[1]), ylab = quote(u[2]),
      zlab = list(deparse(substitute(FUN))[1], rot = 90), \dots)
\S4method{wireframe2}{mvdc}(x, FUN, n.grid = 26, xlim, ylim, zlim = NULL,
      xlab = quote(x[1]), ylab = quote(x[2]),
      zlab = list(deparse(substitute(FUN))[1], rot = 90), \dots)
}
\arguments{
  \item{x}{a \code{"\linkS4class{matrix}"},
    \code{"\linkS4class{data.frame}"}, \code{"\linkS4class{Copula}"} or
    a \code{"\linkS4class{mvdc}"} object.}
  \item{xlim, ylim, zlim}{the x-, y- and z-axis limits.}
  \item{xlab, ylab, zlab}{the x-, y- and z-axis labels.}
  \item{alpha.regions}{see \code{\link{wireframe}()}.}
  \item{scales}{a \code{\link{list}} determining how the axes are drawn;
    see \code{\link{wireframe}()}.}
  \item{par.settings}{See \code{\link{wireframe}()}.}
  \item{FUN}{the \code{\link{function}} to be plotted; typically
    \code{\link{dCopula}} or \code{\link{pCopula}}.}
  \item{n.grid}{the number of grid points used in each dimension.  This
    can be a vector of length two, giving the number of grid points used
    in x- and y-direction, respectively; the function \code{FUN} will be
    evaluated on the corresponding (x,y)-grid.}
  \item{delta}{a small number in \eqn{[0, \frac 1 2)}{[0, 1/2)}
    influencing the evaluation boundaries. The x- and y- vectors will
    have the range \code{[0+delta, 1-delta]}, the default being \code{[0,1]}.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{wireframe}()}, i.e., typically
    \code{\link{panel.3dwire}}, from the \CRANpkg{lattice} package,
    such as \code{shade}, \code{drape}, \code{aspect}, etc.}
}
\section{Methods}{
  Wireframe plots for objects of class \code{"\linkS4class{matrix}"} ,
  \code{"\linkS4class{data.frame}"}, \code{"\linkS4class{Copula}"} or
  \code{"\linkS4class{mvdc}"}.
}
\value{
  An object of class \dQuote{trellis} as returned by \code{\link[lattice]{wireframe}()}.
}
%\author{Marius Hofert}
\seealso{
  The \code{\link{persp-methods}} for drawing perspective plots via base graphics.
}
\examples{
## For 'matrix' objects
## The Frechet--Hoeffding bounds W and M
n.grid <- 26
u <- seq(0, 1, length.out = n.grid)
grid <- expand.grid("u[1]" = u, "u[2]" = u)
W <- function(u) pmax(0, rowSums(u)-1) # lower bound W
M <- function(u) apply(u, 1, min) # upper bound M
x.W <- cbind(grid, "W(u[1],u[2])" = W(grid)) # evaluate W on 'grid'
x.M <- cbind(grid, "M(u[1],u[2])" = M(grid)) # evaluate M on 'grid'
wireframe2(x.W)
wireframe2(x.W, shade = TRUE) # plot of W
wireframe2(x.M, drape = TRUE) # plot of M

## For 'Copula' objects
cop <- frankCopula(-4)
wireframe2(cop, pCopula) # the copula
wireframe2(cop, pCopula, shade = TRUE) # ditto, "shaded"
wireframe2(cop, pCopula, shade = TRUE, col = "gray60") # ditto, "shaded"+grid
wireframe2(cop, pCopula, drape = TRUE, xlab = quote(x[1])) # adjusting an axis label
wireframe2(cop, dCopula, delta=0.01) # the density
wireframe2(cop, dCopula) # => the density is set to 0 on the margins


## For 'mvdc' objects
mvNN <- mvdc(gumbelCopula(3), c("norm", "norm"),
             list(list(mean = 0, sd = 1), list(mean = 1)))
wireframe2(mvNN, dMvdc, xlim=c(-2, 2), ylim=c(-1, 3))
}
\keyword{methods}
\keyword{hplot}
