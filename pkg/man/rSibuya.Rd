\name{rSibuya}
\alias{rSibuya}
\alias{rSibuyaR}
\title{Sampling Sibuya Distributions}
\description{
  Generating random variates from a Sibuya(\eqn{\alpha}{alpha})
  distribution with probability mass function
  \deqn{p_k={\alpha\choose k}(-1)^{k-1},\ k\in\mathbf{N},
  }{p_k = choose(alpha,k)*(-1)^(k-1), k in IN,}
  where \eqn{\alpha\in(0,1]}{alpha in (0,1]}.  The implemented algorithm
  is the one given in Hofert (2011), Proposition 3.2.
}
\usage{
rSibuya(n, alpha)
}
\arguments{
  \item{n}{sample size, i.e., length of the resulting vector of random
    variates.}
  \item{alpha}{parameter in \eqn{(0,1]}.}
}
\value{
  A vector of positive \code{\link{integer}}s of length \code{n} containing the
  generated random variates.
}
\details{
  For documentation and didactical purposes, \code{rSibuyaR} is a pure-\R
  implementation of \code{rSibuya}. However, \code{rSibuyaR} is not as fast as
  \code{rSibuya} (the latter being implemented in C).
}
\author{Marius Hofert, Martin Maechler}
\references{
  Hofert, M. (2011)
  Efficiently sampling nested Archimedean copulas,
  \emph{Computational Statistics & Data Analysis} \bold{55}, 57--70.

  Further, those of the Archimedean families, e.g., \code{\link{copGumbel}}.
}
\examples{
## Sample n random variates from a Sibuya(alpha) distribution and plot a
## histogram
n <- 1000
alpha <- .4
X <- rSibuya(n, alpha)
hist(log(X), prob = TRUE); lines(density(log(X)), col = 2, lwd = 2)
}
\keyword{distribution}
