\name{dDiag}
\alias{dDiag}
\title{Density of Diagonal of an Archimedean Copula}
\description{
  Evaluate the density of the diagonal of a d-dimensional Archimedean copula.
  (Note that the diagonal of a copula is a cumulative distribution function).
}
\usage{
dDiag(u, cop, d, log=FALSE)
}
\arguments{
  \item{u}{numeric vector of evaluation points}
  \item{cop}{an Archimedean acopula object of class
    \code{"\linkS4class{acopula}"}.}
  \item{d}{dimension, an integer \eqn{\ge 1}{>= 1}.}
  \item{log}{logical indicating if \code{\link{log}(f(u))} should be
    returned instead of \code{f(u)}.}
}
\value{numeric vector of density values at \code{u}.}
%\author{Marius Hofert and MM}
% \references{ - Marius?
% }
\seealso{
  \code{\linkS4class{acopula}} class; \code{\link{dnacopula}},
  \code{\link{pnacopula}}.
}
\examples{
th. <- c(0.1, 0.2, 0.5, 0.8, 1.4, 2., 5.)
curve(dDiag(x, cop=setTheta(copClayton, th.[1]), d=3), 0, 1, n=1000,
      ylab = "dDiag(x, *)", main = "diagonal densities of Clayton")
abline(h=0, lty=3)
for(j in 2:length(th.))
  curve(dDiag(x, cop=setTheta(copClayton, th.[j]), d=3), add=TRUE, col=j, n=1000)
legend("topleft", do.call(expression, lapply(th., function(th)
                                 substitute(theta == TH, list(TH=th)))),
       lty = 1, col=seq_along(th.), bty="n")
}
\keyword{distribution}
