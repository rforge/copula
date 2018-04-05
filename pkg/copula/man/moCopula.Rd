\name{moCopula}
\title{The Marshall-Olkin Copula}
\alias{moCopula}
\description{
  Computes Marshall-Olkin copulas in the bivariate case.
}
\usage{
moCopula(param = NA_real_, dim = 2L)
}
\arguments{
  \item{param}{\code{\link{numeric}} vector of length two specifying the
    copula parameters (in \eqn{[0,1]}).}
  \item{dim}{the dimension of the copula.}
}
\value{
  \code{moCopula()} is the constructor for objects of class
  \code{\linkS4class{moCopula}}.
}
\details{
  Note that Marshall-Olkin copulas are only implemented for \code{dim = 2L}.
}
%\author{Marius Hofert}
\examples{
alpha <- c(0.2, 0.8)
MO <- moCopula(alpha)
tau(MO)
lambda(MO)
wireframe2(MO, FUN = pCopula) # if you look carefully, you can see the kink
set.seed(271)
plot(rCopula(1000, MO))
}
\keyword{multivariate}
