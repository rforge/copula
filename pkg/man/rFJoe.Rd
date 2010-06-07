\name{rFJoe}
\alias{rFFrank}
\alias{rFFrankR}
\alias{rFJoe}
\alias{rFJoeR}
\alias{rlog}
\alias{rlogR}
\title{Sample from Copula Related Univariate Distributions}
\description{
  \code{rFFrank}:
  Generate a vector of variates V ~ F with Laplace-Stieltjes transform
  (1-(1-exp(-t)*(1-e^(-theta1)))^alpha)/(1-e^(-theta0)).

  \code{rFJoe}:
  Sample V ~ F with F(n) = 1-1/(n*B(n,1-alpha)), n in IN, with
  Laplace-Stieltjes transform \eqn{1-(1-e^{-t})^\alpha}{1-(1-exp(-t))^alpha}
  via the   algorithm of Hofert (2010).

  \code{rlog}:
  Random number generator for a Log(p) distribution;
  using the algorithm \code{"LK"} of Kemp (1981).
}

\usage{
rFFrank(n, theta0, theta1)
  rFJoe(n, alpha)
   rlog(n, p)
}
\arguments{
  \item{n}{sample size, i.e., length of resulting vector of random variates.}

  \item{theta0}{parameter \eqn{\theta_0 \in (0,\infty)}{theta0 in (0,Inf)}.}
  \item{theta1}{parameter \eqn{\theta_1 \in (\theta_0,\infty)}{theta1 in (theta0,Inf)}.}

  \item{alpha}{parameter}

  \item{p}{parameter in \eqn{(0,1)}.}
}
\details{
  For documentation (and didactical) purposes, all three functions are
  also provided in pure-\R implementation, \code{rFFrankR()},
  \code{rFJoeR()}, and \code{rlogR()}.

...
}
\value{
  a vector \code{n} values; note that these are integers from
  \eqn{1,2,\ldots}{1,2,....}.
}
\references{
  Kemp  ........... (1981).

 ...
}
\seealso{
  \code{\link{rnacopula}},
  \code{\link{rexp}}
}
\examples{
## Frank's F :
## ... %FIXME

## Joe's F :
set.seed(17); x1 <- rFJoe (1000, .4)
set.seed(17); x2 <- rFJoeR(1000, .4)
stopifnot(all.equal(x1,x2))

hist(log(x1), prob=TRUE); lines(density(log(x1)), col=2, lwd=2)

## Extreme Transformation:
summary(llx <- log(1+ log(x1)))
hist(llx, prob=TRUE); lines(density(llx), col=2, lwd=2)

## Log(p) distribution:
## ... %FIXME
}
\keyword{distribution}
