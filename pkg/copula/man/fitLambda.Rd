\name{fitLambda}
\alias{fitLambda}
\title{Non-parametric Estimators of the Matrix of Tail-Dependence Coefficients}
\description{
  Computing non-parametric estimators of the (matrix of) tail-dependence coefficients.
}
\usage{
fitLambda(u, p, tail = c("lower", "upper"), verbose = FALSE)
}
\arguments{
  \item{u}{\eqn{n\times d}{n x d}-\code{\link{matrix}} of (pseudo-)observations in
    \eqn{[0,1]^d} for estimating the (matrix of) tail-dependence coefficients.}
  \item{p}{cut-off parameter in \eqn{[0,1]} below (for \code{tail = "lower"}) or above
  (for \code{tail = "upper"}) which the estimation takes place.}
  \item{tail}{a \code{\link{character}} string specifying whether the
    lower or upper (matrix of) tail-dependence coefficients shall be
    computed.}
  \item{verbose}{a \code{\link{logical}} indicating whether a progress
    bar is displayed.}
}
\value{
  For \code{d = 2} the coefficient of tail dependence. For \code{d > 2},
  the matrix of pairwise coefficients of tail dependence.
}
\details{
  Currently the method of Jaworksi et al. (2009, p. 231) is used; see
  also Schmid and Schmidt (2007).

  As seen in the examples, be careful using such estimators, they might
  not perform too well (depending on \code{p} and in general).
}
%\author{Marius Hofert}
\references{
  Jaworski, P., Durante, F., \enc{Härdle}{Haerdle}, W. K., Rychlik, T. (2010).
  \emph{Copula Theory and Its Applications}
  Springer, Lecture Notes in Statistics -- Proceedings.
}
\examples{
n <- 10000 # sample size
p <- 0.01 # cut-off

## Bivariate case
d <- 2
cop <- claytonCopula(2, dim = d)
set.seed(271)
U <- rCopula(n, copula = cop) # generate observations (unrealistic)
(lam.true <- lambda(cop)) # true tail-dependence coefficients lambda
(lam.C <- c(lower = fitLambda(U, p = p, tail = "lower"),
            upper = fitLambda(U, p = p, tail = "upper"))) # estimate lambdas
## => pretty good
U. <- pobs(U) # pseudo-observations (realistic)
(lam.C. <- c(lower = fitLambda(U., p = p, tail = "lower"),
             upper = fitLambda(U., p = p, tail = "upper"))) # estimate lambdas
## => The pseudo-observations do have an effect...

## Higher-dimensional case
d <- 5
cop <- claytonCopula(2, dim = d)
set.seed(271)
U <- rCopula(n, copula = cop) # generate observations (unrealistic)
(lam.true <- lambda(cop)) # true tail-dependence coefficients lambda
(Lam.C <- list(lower = fitLambda(U, p = p, tail = "lower"),
               upper = fitLambda(U, p = p, tail = "upper"))) # estimate Lambdas
## => Not too good
U. <- pobs(U) # pseudo-observations (realistic)
(Lam.C. <- list(lower = fitLambda(U., p = p, tail = "lower"),
                upper = fitLambda(U., p = p, tail = "upper"))) # estimate Lambdas
## => Performance not too great here in either case
}
\keyword{nonparametric}
\keyword{multivariate}
