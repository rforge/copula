\name{Sibuya}
\alias{rSibuya}
\alias{pSibuya}
\alias{dSibuya}
\alias{rSibuyaR}
\alias{dsumSibuya}
\title{Sibuya Distribution - Sampling and Probabilities}
\description{
  The Sibuya(\eqn{\alpha}{alpha}) distribution can be defined by
  its Laplace transform
  \deqn{1-(1-\exp(-t))^\alpha}{1-(1-exp(-t))^alpha} or its
  probability mass function
  \deqn{p_k={\alpha\choose k}(-1)^{k-1},\ k\in\mathbf{N},
  }{p_k = choose(alpha,k)*(-1)^(k-1), k in IN,}
  where \eqn{\alpha\in(0,1]}{alpha in (0,1]}.

  \code{rSibuya()} generates random variates from
  \eqn{\mathrm{Sib}(\alpha)}{Sib(alpha)},
  by the algorithm given in Hofert (2011), Proposition 3.2.

  \code{dsumSibuya} gives the probability mass function of a sum of
  \eqn{n} Sibuya variates, \eqn{S = \sum_{i=1}^n X_i}{S = sum(i=1..n) X[i]},
  where \eqn{X_i \sim \mathrm{Sib}(\alpha)}{X[i] ~ Sib(alpha)}, i.i.d.
}
\usage{
rSibuya(n, alpha)
dSibuya(x, alpha, log = FALSE)
pSibuya(x, alpha, lower.tail = TRUE, log.p = FALSE)

dsumSibuya(x, n, alpha,
	   method=c("log", "direct", "Rmpfr", "exp.log"), log=FALSE)
}
\arguments{
  \item{n}{
    for \code{rSibuya:} sample size, i.e., length of the resulting
    vector of random variates.
    \cr
    for \code{dsumSibuya}: the number \eqn{n} of summands.
  }
  \item{alpha}{parameter in \eqn{(0,1]}.}
  \item{x}{vector of integer values (\dQuote{quantiles}) at which to
    compute the probability mass or cumulative probability.}
  % \item{p}{vector of probabilities.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
  \item{method}{character string specifying which computational method is to be
    applied.
%% FIXME: explain a bit
    % ##'        log:      proper log computation based on lssum
    % ##'        direct:   brute-force evaluation of the sum and its log
    % ##'        exp.log:  similar to method = "log", but without *proper/intelligent* log
  }
}
\value{
  \describe{% uses bold items --> no \code{..} in there
    \item{rSibuya():}{A vector of positive \code{\link{integer}}s of
      length \code{n} containing the generated random variates.}
    \item{dSibuya(), pSibuya():}{a vector of
      probabilities of the same length as \code{x}.}
    \item{dsumSibuya():}{a vector of probabilities, positive iff
      \code{x >= n} and of the same length as \code{x} (or \code{n} if
      that is longer).

      The value,
      \deqn{\sum_{j=1}^n (-1)^{x-j} {n \choose j}{{j\cdot \alpha}\choose x},%
      }{sum(j=1..n) choose(n,j) * choose(j*alpha,x) * (-1)^(x-j),}
      is numerically difficult to compute for small
      \eqn{\alpha}{alpha} values (even for e.g., \eqn{n \ge 10}{n >=
	10}), and needs high-precision arithmetic there, i.e.,
      \code{method = "Rmpfr"}.}
  }
}
\details{
  The Sibuya distribution has \bold{no} finite moments, i.e., specifically
  infinite mean and variance.

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
