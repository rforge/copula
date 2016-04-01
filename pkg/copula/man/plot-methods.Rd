\name{plot-methods}
\alias{plot-methods}
\alias{plot,copula,ANY-method}
\alias{plot,mvdc,ANY-method}
\docType{methods}
\title{Methods `plot' in Package `copula'}
\description{
  Methods \code{plot()} to draw scatter-plot matrices of
  distributions from package \pkg{copula}.
}
\usage{
\S4method{plot}{copula,ANY}(x, n, xlim = 0:1, ylim = 0:1,
      xlab = expression(U[1]), ylab = expression(U[2]), \dots)
\S4method{plot}{mvdc,ANY}(x, n, xlim = NULL, ylim = NULL,
      xlab = expression(X[1]), ylab = expression(X[2]), \dots)
}
\arguments{
  \item{x}{A \code{"\linkS4class{matrix}"},
    \code{"\linkS4class{data.frame}"}, \code{"\linkS4class{copula}"} or
    a \code{"\linkS4class{mvdc}"} object.}
  \item{n}{The sample size.}
  \item{xlim, ylim}{The x- and y-axis limits.}
  \item{xlab, ylab}{The x- and y-axis labels.}
  \item{\dots}{Additional arguments passed to the underlying
    \code{\link{plot}()}.}
}
\value{invisible().
}
\author{Marius Hofert}
\seealso{
  \code{\link{splom2}()} for a scatter-plot matrix based on \code{\link[lattice]{splom}()}.
}
\examples{
## For 'copula' objects
n <- 1000 # sample size
nu <- 3 # degrees of freedom
tau <- 0.5 # Kendall's tau
th <- iTau(tCopula(df = nu), tau) # corresponding parameter
cop <- tCopula(th, df = nu) # define copula object
set.seed(271)
plot(cop, n = n)

## For 'mvdc' objects
mvNN <- mvdc(cop, c("norm", "norm"),
             list(list(mean = 0, sd = 1), list(mean = 1)))
plot(mvNN, n = n)
}
\keyword{hplot}
