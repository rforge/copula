\name{asymCopula}
\title{Asymmetric (Bivariate or Explicit) Copula Constructors}
\alias{asymBivCopula}
\alias{asymExplicitCopula}
\description{
  Creates an asymmetric copula, an \R object either of
  class \code{"\linkS4class{asymBivCopula}"} or
  class \code{"\linkS4class{asymExplicitCopula}"}.
}
\usage{
asymBivCopula     (copula1, copula2, shapes)
asymExplicitCopula(copula1, copula2, shapes)
}
\arguments{
  \item{copula1, copula2}{each an exchangeable
    \code{\linkS4class{copula}}, of the same dimension \eqn{d},
    bivariate (\eqn{d = 2}) in the case of \code{asymBivCopula()}.}
  \item{shapes}{\code{\link{numeric}} vector of length \eqn{d}, with values
    in \eqn{[0,1]}.}
}
%% \details{
%% }
%\author{Jun Yan and Ivan Kojadinovic}
\value{
  a new \code{"\linkS4class{asymCopula}"} object, an
 \code{"\linkS4class{asymBivCopula}"} or
 \code{"\linkS4class{asymExplicitCopula}"}, respectively.
}
%% \references{
%% }
%% \seealso{
%% }
\examples{
## An asymetric Clayton copula object
ac2 <- asymBivCopula(copula1 = indepCopula(),
                     copula2 = claytonCopula(6),
                     shapes = c(0.4, 0.95))
ac2 ## FIXME: too long 'fullname'
contour(ac2, dCopula, nlevels = 20, main = "dCopula(<asymBivCopula>)")

ac3 <- asymExplicitCopula(copula1 = claytonCopula(2, dim = 3),
                          copula2 = claytonCopula(6, dim = 3),
                          shapes = c(0.2, 0.6, 0.9))
ac3 ## FIXME: too short 'fullname': does not mention 'claytonCopula'
U3 <- rCopula(200, ac3)
pairs2(U3, pch = par("pch"))

}
\keyword{distribution}
\keyword{multivariate}
