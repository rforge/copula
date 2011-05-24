\name{g01}
\alias{g01}
\alias{g01trafo}
\title{Transformation and Anderson-Darling Test of Multivariate Standard Uniform Data}
\description{
  \code{g01()} performs the Anderson-Darling test for \eqn{U[0,1]} of
  supposedly \eqn{U[0,1]^d}-distributed data after it is transformed to
  univariate data via \code{g01trafo()}.

  \code{g01trafo()} transforms supposedly \eqn{U[0,1]^d} distributed vectors of
  random variates to \eqn{U[0,1]}-distributed variates (for checking
  uniformity in a one-dimensional setup).
}
\usage{
g01     (u, method=c("normal","log"))
g01trafo(u, method=c("normal","log"))
}
\arguments{
  \item{u}{\eqn{n\times d}{n x d}-matrix of supposedly \eqn{U[0,1]^d}-distributed data (each value in \eqn{[0,1]}) to be transformed or tested, where \eqn{n} denotes the sample size and \eqn{d} the dimension.}
  \item{method}{a \code{\link{character}} string specifying the method to be used to transform the multivariate data to univariate data.  The method has to be one (or a unique abbreviation) of
    \describe{
		\item{\code{"normal"}}{map to a \eqn{U[0,1]} distribution with the chi-square distribution function with \eqn{d} degrees of freedom applied to the sum of the squares of the quantile function of the standard normal distribution;}
		\item{\code{"log"}}{map to a \eqn{U[0,1]} distribution with the Erlang distribution function (a special gamma distribution) with parameter \eqn{d} applied to the sum of the negative logarithms.}
    }
  }
}
\details{
   For a random vector \eqn{\bm{U}\sim U[0,1]^d}{U~U[0,1]^d} the transformation obtained by \code{method="normal"} returns \eqn{\chi_d^2(\sum_{j=1}^d(\Phi^{-1}(U_j))^2)}{pchisq((Phi^{-1}(U_1))^2+...+(Phi^{-1}(U_d))^2, df=d)} and the one obtained by \code{method="log"} returns \eqn{\Gamma_d(\sum_{j=1}^d(-\log U_j))}{pgamma(-log(U_1)-...-log(U_d), shape=d)}, where \eqn{\chi_d^2}{pchisq(.,df=d)} denotes the distribution function of the chi-square distribution with \eqn{d} degrees of freedom, \eqn{\Phi^{-1}}{Phi^{-1}} denotes the quantile function of the standard normal distribution function, and \eqn{\Gamma_d}{pgamma(.,shape=d)} denotes the distribution function of the gamma distribution with shape parameter \eqn{d} (being equal to an Erlang(\eqn{d}) distribution function).

   Note that the default \code{method="normal"} is normally preferred.
}
\value{
  \code{g01()} returns the Anderson-Darling test result, i.e., the value
  of \code{\link[ADGofTest]{ad.test}()} from package \pkg{ADGofTest}.

  \code{g01trafo()} returns a numeric vector of length \eqn{n}
  (=\code{nrow(u)}) containing the univariate transformed values.
}
\references{
  D'Agostino, R. B., Stephens, M. A. (1986),
  \emph{Goodness-of-fit techniques},
  Dekker.
}
\author{Marius Hofert}
\seealso{
  \code{\link{gnacopula}} (where both functions are applied).
}
\examples{
n <- 200
d <- 20
tau <- 0.25
(theta <- copGumbel@tauInv(tau)) # 4/3
(cop <- onacopulaL("Gumbel", list(theta,1:d)))

## sample
set.seed(1)
U <- rnacopula(n, cop)

## apply the transformation of Hering, Hofert (2011) to obtain U[0,1]^d data
U. <- gnacopulatrafo(U, cop=cop)

## map to U[0,1] data via the chi-square distribution
Z <- g01trafo(U., method="normal")
plot(1:length(Z), Z) # looks uniform as it should

## test via the Anderson-Darling test
g01(U., method="normal") # applies g01trafo
}
\keyword{distribution}
\keyword{multivariate}
