\name{R}
\alias{pR}
\alias{qR}
\title{Distribution of the Radial Part of an Archimedean Copula}
\description{
  \code{pR()} computes the distribution function \eqn{F_R} of the radial
  part of an Archimedean copula. The formula is from McNeil and G.
  \enc{Nešlehová}{Neslehova} (2009).

  \code{qR()} computes the quantile function of \eqn{F_R}. The formula
  is from McNeil and G. \enc{Nešlehová}{Neslehova}
}
\usage{
%%  lower.tail, not survival;  log.p ;  log.p auch fuer qR()
pR(x, family, theta, d, survival = FALSE, log.p = FALSE, \dots)
qR(p, family, theta, d, interval, \dots)
}
\arguments{
  \item{x}{numeric vector of nonnegative evaluation points for \eqn{F_R}.}
  \item{p}{numeric vector %, in \eqn{[0,1]},
    of evaluation points of the  quantile function.}
  \item{family}{Archimedean family.}
  \item{theta}{parameter \eqn{theta}.}
  \item{d}{dimension \eqn{d}.}
  \item{survival, log}{\code{\link{logical}} indicating whether the
    \emph{logarithm} of the \emph{survival} function .. cfunction of \eqn{R} is
    computed.}
  \item{interval}{root-search interval.}
  \item{\dots}{additional arguments passed to the procedure for
    computing derivatives.}
}
\value{The distribution function of the radial part evaluated at
  \code{x}, or its inverse, the quantile at \code{p}.}
\author{Marius Hofert}
\references{
  McNeil, A. J., G. \enc{Nešlehová}{Neslehova}, J. (2009).
  Multivariate Archimedean copulas, \eqn{d}-monotone functions and
  \eqn{l_1}-norm symmetric distributions. \emph{The Annals of Statistics}
  \bold{37(5b)}, 3059--3097.
}
\examples{
## setup
family <- "Gumbel"
tau <- 0.5
m <- 256
dmax <- 20
x <- seq(0, 20, length.out=m)

## compute and plot pR() for various d's
y <- vapply(1:dmax, function(d)
            pR(x, family=family, theta=iTau(archmCopula(family), tau), d=d),
            rep(NA_real_, m))
plot(x, y[,1], type="l", ylim=c(0,1),
     xlab = expression(italic(x)),
     ylab = substitute(italic(F[R](x))~~"for d=1:"*dm, list(dm=dmax)))
for(k in 2:dmax) lines(x, y[,k])
}
\keyword{distribution}
