\name{math-fun}
\alias{A..Z}
\alias{sinc}
\title{Sinc, Zolotarev's, and other mathematical utility functions}
\description{
  \code{A..Z(x, a)} computes Zolotarev's function to
  the power \code{1-a}.

  \code{sinc(x)} computes the \dQuote{sinc} function
  \eqn{s(x) = sin(x)/x} for \eqn{x\ne 0}{x != 0} and \eqn{s(0) = 1},
  such that \eqn{s()} is continuous, also at \eqn{x = 0}.
}
\usage{
A..Z(x, alpha, I.alpha = 1 - alpha)

sinc(x)
}
\arguments{
  \item{x}{numeric argument, typically a vector.}
  \item{alpha}{parameter in (0,1].}
  \item{I.alpha}{must be \code{ = 1 - alpha}, maybe more accurately
    when \code{alpha} is very close to 1.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \code{A..Z(x,alpha)} is \eqn{\tilde A_{Z}(x,\alpha)}{A~Z(x,alpha)},
  defined as
  \deqn{\frac{sin(\alpha x)^\alpha * sin((1-\alpha)x)^(1-\alpha)}{sin(x)},}{sin(alpha*x)^alpha * sin((1-alpha)*x)^(1-alpha) / sin(x),}

  where \eqn{\alpha}{alpha} is \code{alpha}.
}
%\author{Martin Maechler}
\references{
  ...

  Luc Devroye (2010) ...
  ... % FIXME
}
\seealso{
  \code{\link{retstable}} internally makes use of these functions.
}
\examples{
curve(sinc, -15,25); abline(h=0,v=0, lty=2)
curve(A..Z(x, 0.25), -4,4,
      main = "Zolotarev's function A(x) ^ 1-alpha")
}
\keyword{math}
