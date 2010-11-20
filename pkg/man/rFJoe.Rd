\name{rFJoe}
\alias{rFFrank}
\alias{rFFrankR}
\alias{rFJoe}
\alias{rFJoeR}
\alias{rlog}
\alias{rlogR}
\title{Sample from Univariate Distributions involved in nested Frank and Joe Copulas}
\description{
  \code{rlog}: Generating random variates of a Log(p) distribution with probability mass function \deqn{p_k=p(1-p)^{k-1},\ k\in\mathbf{N},}{p_k = p(1-p)^(k-1), k in IN,} where \eqn{p\in(0,1)}{p in (0,1)}. The implemented algorithm is the one named \dQuote{LK} in Kemp (1981).

  \code{rFJoe}: Generating random variates from the distribution with probability mass function \deqn{p_k={\alpha\choose k}(-1)^{k-1},\ k\in\mathbf{N},}{p_k = binom{alpha}{k}(-1)^(k-1), k in IN,} where \eqn{\alpha\in(0,1]}{alpha in (0,1]}. This distribution appears on sampling (nested) Joe copulas. Here, \eqn{F_0}{F0} and \eqn{F_{01}}{F01} belong to this class. The implemented algorithm is given in Hofert (2010a).

  \code{rFFrank}: Generating random variates from the distribution with probability mass function \deqn{p_k={\theta_0/\theta_1\choose k}(-1)^{k-1}(1-e^{-\theta_1})^k/(1-e^{-\theta_0}),\ k\in\mathbf{N},}{p_k = binom{theta0/theta1}{k}(-1)^(k-1)(1-e^(-theta1))^k/(1-e^(-theta0)), k in IN,} where \eqn{\theta_0,\theta_1\in(0,\infty)}{theta0, theta1 in (0,Inf)} such that \eqn{\theta_0\le\theta_1}{theta0 <= theta1}. This distribution appears on sampling nested Frank copulas. Here, \eqn{F_{01}}{F01} is a \eqn{V_0}{V0}-fold convolution of random variables from this distribution. The implemented algorithm is given in Hofert (2010a).
}

\usage{
rlog(n, p, Ip = 1 - p)
rFJoe(n, alpha)
rFFrank(n, theta0, theta1)
}
\arguments{
  \item{n}{sample size, i.e., length of the resulting vector of random variates.}
  \item{p}{parameter in \eqn{(0,1)}.}
  \item{Ip}{\eqn{= 1 - p}, possibly more accurately, e.g, when
    \eqn{p\approx 1}{p ~= 1}.}
  \item{alpha}{parameter in \eqn{(0,1]}.}
  \item{theta0, theta1}{parameters \eqn{\theta_0,\theta_1\in (0,\infty)}{theta0, theta1 in (0,Inf)} such that \eqn{\theta_0\le\theta_1}{theta0 <= theta1}.}
}
\value{
  A vector of positive \code{\link{integer}}s of length \code{n} containing the generated random variates.
}
\details{
  For documentation (and didactical) purposes, all three functions are
  also provided in pure-\R implementation, \code{rlogR()},
  \code{rFJoeR()}, and \code{rFFrankR()}.
}
\author{Marius Hofert, Martin Maechler}
\seealso{
  \code{\link{rnacopula}}
}
\references{
  	Kemp, A. W. (1981),
	Efficient Generation of Logarithmically Distributed Pseudo-Random Variables,
	\emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)} \bold{30}, 3, 249--253.

	Hofert, M. (2010a),
	Efficiently sampling nested Archimedean copulas,
	\emph{Computational Statistics & Data Analysis}, in press.

	Hofert, M. (2010b),
	\emph{Sampling Nested Archimedean Copulas with Applications to CDO Pricing},
	Suedwestdeutscher Verlag fuer Hochschulschriften AG & Co. KG.
}
\examples{
## Sample n random variates from a Log(p) distribution and plot a
## histogram
n <- 500
p <- .5
X <- rlog(n, p)
hist(X, prob=TRUE)

## Sample the univariate distribution appearing when sampling nested
## Joe copulas and plot a histogram
alpha <- .4
Y <- rFJoe(n, alpha)
hist(log(Y), prob=TRUE); lines(density(log(Y)), col=2, lwd=2)

## Sample the univariate distribution appearing when sampling nested
## Frank copulas and plot a histogram
theta0 <- 3
theta1 <- 8
Z <- rFFrank(n, theta0, theta1)
hist(log(Z), prob=TRUE); lines(density(log(Z)), col=2, lwd=2)
}
\keyword{distribution}
