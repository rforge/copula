\name{rotCopula}
\alias{rotCopula}

\title{Construction of "Rotated" Copulas}
\description{
  Constructs a rotated copula from an initial copula and a vector of
  logicals indicating which dimension to \dQuote{flip}.
}
\usage{
rotCopula(copula, flip = rep(TRUE, copula@dimension))
}
\arguments{
  \item{copula}{an object of class \code{"\linkS4class{Copula}"}.}
  \item{flip}{a vector of logicals indicating which dimension should be
    \dQuote{flipped}; by default, all the components are rotated, implying that
    the rotated copula is the survival copula.}
}
\value{
  A rotated copula object of class \code{"rotCopula"}.
}
%\details{}
%\author{Ivan Kojadinovic}
%\references{}
\seealso{
  \code{"\linkS4class{Copula}"},
  \code{\link{fitCopula}} for fitting such copulas to data.
}
\examples{
## A two-dimensional example: a "rotated" Clayton copula
rc <- rotCopula(claytonCopula(3), flip = c(TRUE, FALSE))

contour(rc, dCopula, nlevels = 20)
contour(rc, pCopula, nlevels = 20)
rho(rc)
tau(rc)

n <- 1000
u <- rCopula(n, rc)
rho.n <- cor(u[,1], u[,2], method = "spearman")
tau.n <- cor(u[,1], u[,2], method = "kendall")

iRho(rotCopula(claytonCopula(), flip = c(TRUE, FALSE)), rho.n)
iTau(rotCopula(claytonCopula(), flip = c(TRUE, FALSE)), tau.n)

## Fitting
fitCopula(rotCopula(claytonCopula(), flip = c(TRUE, FALSE)), pobs(u),
method = "irho")
fitCopula(rotCopula(claytonCopula(), flip = c(TRUE, FALSE)), pobs(u),
method = "itau")
fitCopula(rotCopula(claytonCopula(), flip = c(TRUE, FALSE)), pobs(u),
method = "mpl")

\donttest{
## Goodness-of-fit testing
gofCopula(rotCopula(claytonCopula(), flip = c(TRUE, FALSE)), u)
gofCopula(rotCopula(claytonCopula(), flip = c(TRUE, FALSE)), u,
          sim = "mult")
}

## A four-dimensional example: a "rotated" Frank copula
rf <- rotCopula(frankCopula(10, dim = 4),
                flip = c(TRUE, FALSE, TRUE, FALSE))

n <- 1000
u <- rCopula(n, rf)
pairs(u)

pCopula(c(0.6,0.7,0.6,0.8), rf)
C.n(matrix(c(0.6,0.7,0.6,0.8), 1, 4), u)

## Fitting: itau and irho should not be used (FIXME?)
fitCopula(rotCopula(frankCopula(dim=4),
                    flip = c(TRUE, FALSE, TRUE, FALSE)), pobs(u))

\dontrun{
## Goodness-of-fit testing
gofCopula(rotCopula(frankCopula(dim=4),
                    flip = c(TRUE, FALSE, TRUE, FALSE)), pobs(u))
gofCopula(rotCopula(frankCopula(dim=4),
                    flip = c(TRUE, FALSE, TRUE, FALSE)), pobs(u),
          sim = "mult")
}
}
\keyword{distribution}
\keyword{multivariate}
