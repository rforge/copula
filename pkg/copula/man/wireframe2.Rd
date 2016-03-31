\name{wireframe2}
\alias{wireframe2}
\title{Wireframe Plot with Nice Defaults}
\description{
  A version of \pkg{lattice}'s \code{\link{wireframe}} function,
  particularly useful for visualizing copulas or their densities.
}
\usage{
wireframe2(x, labels = NULL,
           labels.null.lab = parse(text = c("u[1]", "u[2]", "C(u[1],u[2])")),
           alpha.regions = 0.5, scales = list(arrows = FALSE, col = 1),
           par.settings = standard.theme(color = FALSE), ...)
}
\arguments{
  \item{x}{A numeric matrix or as.matrix(.)able.}
  \item{labels}{The x-, y- and z-axis labels.}
  \item{labels.null.lab}{A character vector of length 3 determining the
    \dQuote{base name} of the variable labels in case \code{x} does not
    have all column names given and \code{labels} is \code{NULL}.}
  \item{alpha.regions}{See \code{\link{wireframe}()}.}
  \item{scales}{See \code{\link{wireframe}()}.}
  \item{par.settings}{See \code{\link{wireframe}()}.}
  \item{\dots}{Further arguments, passed to \code{\link[lattice]{wireframe}()}.}
}
\value{From \code{\link{wireframe}()}, an \R object of class \code{"trellis"}.
}
\author{Marius Hofert.}
\examples{
## Wireframe plots of the Frechet--Hoeffding lower and upper bounds
W <- function(u) pmax(0, rowSums(u)-1) # lower bound W
M <- function(u) apply(u, 1, min) # upper bound M
m <- 20 # grid size is m x m
u <- seq(0, 1, length.out = m) # subdivision points in each dimension
grid <- expand.grid("u[1]" = u, "u[2]" = u) # grid
W.grid <- cbind(grid, "W(u[1],u[2])" = W(grid)) # evaluate W on the grid
M.grid <- cbind(grid, "M(u[1],u[2])" = M(grid)) # evaluate M on the grid
wireframe2(W.grid) # plot of W
wireframe2(M.grid) # plot of M
}
\keyword{hplot}
