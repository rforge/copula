\name{rotCopula}
\title{Construction and Class of "Rotated" Copulas}
\alias{rotCopula}
\alias{rotCopula-class}
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
    \dQuote{flipped}; by default, all the components are flipped, implying that
    the rotated copula is the survival copula.}
}
\value{
  A rotated copula object of class \code{"rotCopula"}.
}
\section{Slots}{
  of a \code{"rotCopula"} object
  \describe{
    \item{\code{copula}:}{Object of class \code{"\linkS4class{copula}"}.}
    \item{\code{flip}:}{\code{\link{logical}} vector of length \eqn{d}
      (the copula dimension) specifying which margins are flipped;
      corresponds to the \code{flip} argument of \code{rotCopula()}. }
    \item{\code{dimension}:}{the copula dimension \eqn{d}, an \code{\link{integer}}.}
    \item{\code{parameters}:}{\code{\link{numeric}} vector specifying
      the parameters.}
    \item{\code{param.lowbnd}, and \code{param.upbnd}:}{\code{numeric}
      vector of the same length as \code{parameters}, specifying
      (component wise) bounds for each of the parameters.}
    \item{\code{param.names}:}{\code{\link{character}} vector (of same
      length as \code{parameters}) with parameter names.}
    \item{\code{fullname}:}{a \code{character} string describing the
      rotated copula.}
  }
}
%\details{}
%\author{Ivan Kojadinovic}
%\references{}
\seealso{
  \code{\link{fitCopula}} for fitting such copulas to data.
}
\examples{
## A two-dimensional example: a "rotated" Clayton copula
rc <- rotCopula(claytonCopula(3), flip = c(TRUE, FALSE))

contour(rc, dCopula, nlevels = 20, main = "dCopula(<rotCopula>)")
contour(rc, pCopula, nlevels = 20, main = "pCopula(<rotCopula>)")
rho(rc)
tau(rc) # -0.6

n <- 1000
u <- rCopula(n, rc)
rho.n <- cor(u[,1], u[,2], method = "spearman")
tau.n <- cor(u[,1], u[,2], method = "kendall")

## "Calibration"
rc. <- rotCopula(claytonCopula(), flip = c(TRUE, FALSE))
iRho(rc., rho.n)
iTau(rc., tau.n)

## Fitting
fitCopula(rc., pobs(u), method = "irho")
fitCopula(rc., pobs(u), method = "itau")
fitCopula(rc., pobs(u), method = "mpl")

## Goodness-of-fit testing -- the first, parametric bootstrap, is *really* slow
\dontrun{gofCopula(rc., u)}
\donttest{gofCopula(rc., u, sim = "mult")}

## A four-dimensional example: a "rotated" Frank copula
rf <- rotCopula(frankCopula(10, dim = 4),
                flip = c(TRUE, FALSE, TRUE, FALSE))

n <- 1000
u <- rCopula(n, rf)
pairs(u)

pCopula(c(0.6,0.7,0.6,0.8), rf)
C.n(cbind(0.6,0.7,0.6,0.8),  u)

## Fitting: itau and irho should not be used (FIXME?)
(rf. <- rotCopula(frankCopula(dim=4),
                  flip = c(TRUE, FALSE, TRUE, FALSE)))
fitCopula(rf., pobs(u))

## Goodness-of-fit testing (first ("PB") is really slow, see above):
\dontrun{gofCopula(rf., pops(u))}
\donttest{gofCopula(rf., pobs(u), sim = "mult") # takes 3.7 sec [lynne, 2015]}
}
\keyword{distribution}
\keyword{multivariate}
