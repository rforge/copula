\name{rotCopula}
\alias{rotCopula}

\title{Construction of Rotated Copulas}
\description{
  Constructs a rotated copula from an initial copula and a mask of
  logicals.
}
\usage{
rotCopula <- function(copula, mask = rep(TRUE, copula@dimension))
}
\arguments{
  \item{copula}{an object of class \code{"\linkS4class{Copula}"}.}
  \item{mask}{a vector of logicals indicating which component should be
    rotated; by default, all the components are rotated, implying that
    the rotated copula is the survival copula.}
}
\value{
  A Rotated copula object of class \code{"rotCopula"}.
}
%\details{}
%\references{}
\seealso{
  \code{\linkS4class{Copula}},
  \code{\link{fitCopula}} for fitting such copulas to data.
}
\examples{
\donttest{
## A two-dimensional example
\dontrun{
rotated.clayton <- rotCopula(claytonCopula(3), mask = c(T, F))

contour(rotated.clayton, dCopula, nlevels = 20)
contour(rotated.clayton, pCopula, nlevels = 20)
rho(rotated.clayton)
tau(rotated.clayton)

n <- 1000
u <- rCopula(n, rotated.clayton)
rho.n <- cor(u[,1], u[,2], method = "spearman")
tau.n <- cor(u[,1], u[,2], method = "kendall")

iRho(rotCopula(claytonCopula(), mask = c(T, F)), rho.n)
iTau(rotCopula(claytonCopula(), mask = c(T, F)), tau.n)

fitCopula(rotCopula(claytonCopula(), mask = c(T, F)), pobs(u),
method = "itau")
fitCopula(rotCopula(claytonCopula(), mask = c(T, F)), pobs(u),
method = "irho")
fitCopula(rotCopula(claytonCopula(), mask = c(T, F)), pobs(u),
method = "mpl")

gofCopula(rotCopula(claytonCopula(), mask = c(T, F)), u)
gofCopula(rotCopula(claytonCopula(), mask = c(T, F)), u, sim = "mult")

## A four-dimensional example
rotated.frank <- rotCopula(frankCopula(10, dim = 4), mask = c(T, F, T, F))
n <- 1000
u <- rCopula(n, rotated.frank)
pairs(u)

pCopula(c(0.6,0.7,0.6,0.8), rotated.frank)
C.n(matrix(c(0.6,0.7,0.6,0.8),1,4), u)


fitCopula(rotCopula(frankCopula(dim=4), mask = c(T, F, T, F)), pobs(u))
## itau and irho should not be used
gofCopula(rotCopula(frankCopula(dim=4), mask = c(T, F, T, F)), pobs(u))
gofCopula(rotCopula(frankCopula(dim=4), mask = c(T, F, T, F)), pobs(u), sim = "mult")
} # dontrun
}
}
\keyword{distribution}
\keyword{multivariate}
