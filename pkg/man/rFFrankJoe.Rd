\name{rFFrankJoe}
\alias{rFFrank}
\alias{rFJoe}
\title{Sampling Distribution F for Frank and Joe}
\description{
  Generate a vector of variates \eqn{V \sim F}{V ~ F} with
  Laplace-Stieltjes transform
  \deqn{(1-(1-exp(-t)*(1-e^(-theta1)))^alpha)/(1-e^(-theta0)),}
  for Frank, or
  \deqn{1-(1-exp(-t))^alpha,}
  for Joe, respectively.
}
\usage{
rFFrank(n, theta0, theta1, rej)
rFJoe(n, alpha)
}
\arguments{
  \item{n}{number of variates from F}
  \item{theta0}{parameter theta0 in (0,infinity)}
  \item{theta1}{parameter theta1 in [theta0, infinity)}
  \item{rej}{method switch: if theta_0 > rej a rejection
    from Joe (Sibuya distribution) is applied (otherwise, a logarithmic
    envelope is used).}
  \item{alpha}{for \code{rFJoe()}, the parameter alpha is \eqn{\alpha=
      \theta_0/\theta_1}{alpha = theta0/theta1}, must lie in \eqn{(0,1]}.}
}
\value{
  numeric vector of of random variates \eqn{V}, of length \code{n}.
}
% \details{
% }
\author{Marius Hofert}
\seealso{
  \code{\link{rF01Frank}} and \code{\link{rF01Joe}};% in ./rF01FrankJoe.Rd
  \code{\link{rnacopula}}.
}
\examples{
## FIXME
}
\keyword{distribution}
