\name{gnacopula}
\alias{gnacopula}
\alias{gnacopulatrafo}
\title{Goodness of Fit Testing for Nested Archimedean Copulas}
\description{
  \code{gnacopula()} conducts a goodness-of-fit test for the given H0
  copula \code{cop} based on the data \code{u}.

  \code{gnacopulatrafo()} transforms vectors of random variates following
  the given nested Archimedean copula (with specified parameters) to
  \eqn{U[0,1]^d} vectors of random variates.
}
\usage{
gnacopula(u, cop, n.bootstrap = 0,
          estimation.method = eval(formals(enacopula)$method),
          include.K = TRUE, n.MC = 0, method = c("normal", "log"),
          verbose = TRUE, \dots)

gnacopulatrafo(u, cop, include.K = TRUE, n.MC = 0)
}
\arguments{
  \item{u}{(pseudo) data matrix, of values in \eqn{[0,1]}.}
  \item{cop}{\code{"\linkS4class{outer_nacopula}"}, for
    \code{gnacopula()} with specified \eqn{H_0}{H0} parameters.}
  \item{n.bootstrap}{number of bootstrap replications.}
  \item{estimation.method}{estimation method, see \code{\link{enacopula}}.}
  \item{include.K}{logical indicating whether the last component,
    \eqn{K()}, is also used or not.}
  \item{n.MC}{parameter \code{n.MC} for \code{\link{K}()}.}
  \item{method}{for \code{\link{g01}}, see there.}
  \item{\dots}{additional arguments to \code{\link{enacopula}}.}
  \item{verbose}{if TRUE, the progress of the bootstrap is displayed.}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  \code{gnacopula()} returns an \R object of class \code{"htest"}.
  In the case of \code{n.bootstrap = 0}, simply
  \code{\link{g01}(gnacopulatrafo(....))}, otherwise the result also
  contains the bootstrap results, and is a list with components
 \item{comp1 }{Description of 'comp1'}
 \item{comp2 }{Description of 'comp2'}

  FIXME

  \code{gnacopulatrafo()} returns a numeric vector of length \code{nrow(u)}.
}
\references{
  FIXME %% ~put references to the literature/web site here ~
}
\author{Marius Hofert and Martin Maechler}
\seealso{
  \code{\link{g01}()} and its transformation.
}
\examples{
 ## Note that
 ##    demo(estimation.gof)
 ##    -------------------- is a whole small scale simulation, using gnacopula()

%%%%%%  FIXME ......................
}
\keyword{htest}
\keyword{distribution}
\keyword{multivariate}

