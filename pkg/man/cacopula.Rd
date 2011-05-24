\name{cacopula}
\alias{cacopula}
\title{Conditional Copula Function}
\description{
  Compute the conditional copula function \eqn{C(v|u)} of
  \eqn{v} given \eqn{u}.
}
\usage{
cacopula(v, u, family, theta, log=FALSE)
}
\arguments{
  \item{v}{parameter v}
  \item{u}{parameter u}
  \item{family}{Archimedean family (name or object)}
  \item{theta}{parameter theta}
  \item{log}{logical indicating if \eqn{\log(C(v|u))} should be retruen is returned.}
}
\value{
  ... FIXME
}
\author{Marius Hofert}
\note{
  for some families, this function makes sense for \eqn{u = 0} or
  \eqn{v = 0}, since the corresponding limits can be computed; but not for all.
}
\seealso{
  \code{\link{acopula-families}}.
}
\examples{
theta <- 2
C2 <- onacopula("Clayton", C(theta, 1:2))
C2@copula # is an "acopula" with specific parameter theta
tt <- (1:15)/16
cacopula(v = tt, u = rev(tt),     "Clayton", theta = 2)
cacopula(v = tt, u = rep(0.5,15), "Clayton", theta = 2)
## and with a wrong theta : % MM to MH: does this make *ANY* sense?
cacopula(v = tt, u = rep(0.5,15), "Clayton", theta = 5)
}
\keyword{distribution}
