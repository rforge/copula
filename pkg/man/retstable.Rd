\name{retstable}
\alias{retstable}
\alias{retstableR}
\title{Generate Exponentially Tilted Stable Random Variates St}
\description{Generate a sample vector of exponentially tilted stable
  random variates \eqn{St ~ \tilde{S}(\alpha, 1,
    (cos(\alpha*pi/2) *V_0)^{1/\alpha},
    V_0*I_{\alpha = 1}, I_{\alpha != 1}; 1)}

  with Laplace-Stieltjes transform
  \eqn{exp(-V_0((1+t)^\alpha-1))}{exp(-V0((1+t)^alpha -1))};

  see Nolan's book for the parametrization.
}
\usage{
retstableR(alpha, V0, h = 1)

retstable(alpha, V0, h = 1, method = NULL)
}
\arguments{
  \item{alpha}{parameter in (0,1]}
  \item{V0}{vector of random variates}
  \item{h}{non-negative real number}
  \item{method}{a character string denoting the method to use, currently
    either \code{"MH"} (Marius Hofert's) or \code{"LD"} (Luc Devroye's).
    By default, when \code{NULL}, a smart choice is made to use the
    faster method.
  }
}
\details{
  \code{retstableR} is version of \code{"MH"}, in a pure \R
  implementation.

  ...
}
\value{
  vector of variates St
}
\references{
  Nolan, ...
}
\seealso{
 \code{\link{rstable1}}, ... % FIXME
}
\examples{
set.seed(2); rETS <- retstable (.2, V0=c(1:200))
# in order to be comparable to *R(), specify method:
set.seed(2); rETS. <- retstable (.2, V0=c(1:200), method = "MH")
## the R version takes  *much* longer:
set.seed(2); rETS2 <- retstableR(.2, V0=c(1:200))
stopifnot(all.equal(rETS., rETS2))
}
\keyword{distribution}
