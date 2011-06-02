\name{tauAMH}
\alias{tauAMH}
\title{Ali-Mikhail-Haq ("AMH")'s  Kendall's Tau}
\description{
  Compute Kendall's Tau of an Ali-Mikhail-Haq ("AMH") Archimedean copula
  with parameter \code{theta}.  While that's analytically given
  explicitly, as
  \deqn{\frac{1 - 2((1-\theta)(1-\theta)\cdot\log(1-\theta) + \theta)}{3\theta^2},}{
    1 - 2((1-t)(1-t)log(1-t) + t)/(3*t^2),}
  for \code{th}\eqn{=\theta}{=t};
  numerically, care has to be taken when \eqn{\theta \to 0}, avoiding accuracy
  loss already, e.g., for \eqn{\theta} as large as \code{theta = 0.001}.
}
\usage{
tauAMH(th)
}
\arguments{
  \item{th}{numeric vector with values in \eqn{[0,1]}.}
}
\value{
  a vector of the same length as \code{th},
  with values of \code{2*((1-th)*(1-th)*log(1-th) + th) / (3*th^2)},
  numerically accurately, to at least around 12 decimal digits.
}
\author{Martin Maechler}
\seealso{
  \code{\link{acopula-families}}, and their class definition,
  \code{"\linkS4class{acopula}"}.
}
\examples{
tauAMH(c(0, 2^-40, 2^-20))
curve(tauAMH, 0, 1)
curve(tauAMH, 1e-12, 1, log="xy") # linear, tau ~= 2/9*theta in the limit
}
\keyword{distribution}
