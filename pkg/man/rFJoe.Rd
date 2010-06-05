\name{rFJoe}
\alias{rFJoe}
\alias{rFJoeR}
\title{Sample Joe' Copula F}
\description{Sample V ~ F with F(n) = 1-1/(n*B(n,1-alpha)), n in IN, with
  Laplace-Stieltjes transform \eqn{1-(1-e^{-t})^\alpha}{1-(1-exp(-t))^alpha}
  via the   algorithm of Hofert (2010).
}

\usage{
rFJoe(n, alpha)
}
\arguments{
  \item{n}{sample size}
  \item{alpha}{parameter}
}
\details{
...
}
\value{
  a vector \code{n} values; note that these are integers from
  \eqn{1,2,\ldots}{1,2,....}.
}
\references{
 ...
}
\seealso{
  \code{\link{rnacopula}},
  \code{\link{rexp}}
}
\examples{
set.seed(17); x1 <- rFJoe (1000, .4)
set.seed(17); x2 <- rFJoeR(1000, .4)
stopifnot(all.equal(x1,x2))

hist(log(x1), prob=TRUE); lines(density(log(x1)), col=2, lwd=2)

## Extreme Transformation:
summary(llx <- log(1+ log(x1)))
hist(llx, prob=TRUE); lines(density(llx), col=2, lwd=2)
}
\keyword{distribution}
