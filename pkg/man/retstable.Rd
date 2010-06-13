\name{retstable}
\alias{retstable}
\alias{retstableR}
\title{Sampling Exponentially Tilted Stable Random Variates}
\description{Generating random variates of an exponentially tilted stable 
distribution of the form \deqn{\tilde{S}(\alpha, 1,
    (\cos(\alpha\pi/2)V_0)^{1/\alpha},
    V_0\mathbf{1}_{\{\alpha=1\}}, h\mathbf{1}_{\{\alpha\neq1\}}; 1)}{tS(alpha, 1, (cos(alpha*pi/2)V0)^(1/alpha), V0*1_(alpha==1), h*1_(alpha!=1))}
  with parameters \eqn{\alpha\in(0,1]}{alpha in (0,1]}, \eqn{V_0\in(0,\infty)}{V0 in (0,Inf)}, and \eqn{h\in[0,\infty)}{h in [0,Inf)} and corresponding Laplace-Stieltjes transform
  \deqn{\exp(-V_0((h+t)^\alpha-h^\alpha)),\ t\in[0,\infty],}{exp(-V0((h+t)^alpha-h^alpha)), t in [0,Inf],} see the references for more details about this distribution.
}
\usage{
retstable(alpha, V0, h = 1, method = NULL)
retstableR(alpha, V0, h = 1)
}
\arguments{
  \item{alpha}{parameter in \eqn{(0,1]}.}
  \item{V0}{vector of values in \eqn{(0,\infty)}{(0,Inf)} (e.g., when sampling nested Clayton copulas, these are random variates from \eqn{F_0}{F0}).}
  \item{h}{parameter in \eqn{[0,\infty)}{[0,Inf)}.}
  \item{method}{a character string denoting the method to use, currently
    either \code{"MH"} (Marius Hofert's algorithm) or \code{"LD"} (Luc Devroye's algorithm). By default, when \code{NULL}, a smart choice is made to use the
    faster method depending on the specific values of \eqn{V_0}{V0}.}
  }
}
\details{
  \code{retstableR} is version of \code{"MH"}, in a pure \R
  implementation, however not as fast as \code{retstable} and therefore not recommended in simulations when run time matters.
}
\value{
  A vector of variates from \eqn{\tilde{S}(\alpha, 1,
    (\cos(\alpha\pi/2)V_0)^{1/\alpha},
    V_0\mathbf{1}_{\{\alpha=1\}}, h\mathbf{1}_{\{\alpha\neq1\}}; 1)}{tS(alpha, 1, (cos(alpha*pi/2)V0)^(1/alpha), V0*1_(alpha==1), h*1_(alpha!=1))}.
}
\author{Marius Hofert, Martin Maechler}
\seealso{
 \code{\link{rstable1}} for sampling stable distributions.
}
\references{
	Hofert, M. (2010a),
	Efficiently sampling nested Archimedean copulas,
	\emph{Computational Statistics & Data Analysis}, in press.
	
	Hofert, M. (2010b),
	\emph{Sampling Nested Archimedean Copulas with Applications to CDO Pricing},
	Suedwestdeutscher Verlag fuer Hochschulschriften AG & Co. KG.
}
\examples{
## Draw random variates from an exponentially tilted stable distribution
## with given alpha, V0, and h = 1
alpha <- .2
V0 <- rgamma(200, 1)
rETS <- retstable(alpha, V0)

## Plot the random variates
plot(rETS)
}
\keyword{distribution}
