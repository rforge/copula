\name{fitLambda}
\alias{fitLambda}
\title{Non-parametric Estimators of the Matrix of Tail-Dependence Coefficients}
\description{
  Computing non-parametric estimators of the (matrix of) tail-dependence coefficients.
}
\usage{
fitLambda(u, method = c("Schmid.Schmidt", "Schmidt.Stadtmueller", "t"),
          p = 1/sqrt(nrow(u)), lower.tail = TRUE, verbose = FALSE, ...)
}
\arguments{
  \item{u}{\eqn{n\times d}{n x d}-\code{\link{matrix}} of (pseudo-)observations in
    \eqn{[0,1]^d} for estimating the (matrix of) tail-dependence coefficients.}
  \item{method}{the method with which the tail-dependence coefficients
    are computed:
    \describe{
      \item{\code{method = "Schmid.Schmidt"}:}{nonparametric estimator
	of Schmid and Schmidt (2007) (see also Jaworksi et al. (2009, p. 231))
        computed for all pairs (pairwise conditional Spearman's rho
	for \code{u <= p}).}
      \item{\code{method = "Schmidt.Stadtmueller"}:}{nonparametric
        estimator of Schmidt and Stadtmueller (2006)
        computed for all pairs (pairwise empirical tail copula
	for \code{u <= p}).}
      \item{\code{method = "t"}:}{fits pairwise \eqn{t} copulas and
	returns the implied tail-dependence coefficient.}
    }
  }
  \item{p}{(small) cut-off parameter in \eqn{[0,1]} below (for
    \code{tail = "lower"}) or above (for \code{tail = "upper"}) which
    the estimation takes place.}
  \item{lower.tail}{\code{\link{logical}} indicating whether the lower
    (the default) or upper tail-dependence coefficient is computed; in
    case of the latter, the lower tail dependence coefficient of
    the flipped data \code{1-u} is computed.}
  \item{verbose}{a \code{\link{logical}} indicating whether a progress
    bar is displayed.}
  \item{\dots}{additional arguments passed to the underlying functions
    (at the moment only to \code{\link{optimize}()} in case \code{method
    = "t"}).}
}
\value{
  The matrix of pairwise coefficients of tail dependence; for
  \code{method = "t"} a \code{\link{list}} with components
  \code{Lambda},
  the matrix of pairwise estimated correlation coefficients \code{P}
  and the matrix of pairwise estimated degrees of freedom \code{Nu}.
}
\details{
  As seen in the examples, be careful using nonparametric estimators,
  they might not perform too well (depending on \code{p} and in
  general). After all, the notion of tail dependence is a limit,
  finite sample sizes may not be able to capture limits well.
}
%\author{Marius Hofert}
\seealso{\code{\link{lambda}()} computes the true (lower and upper) tail
  coefficients for a given copula.
}
\references{
  Jaworski, P., Durante, F., \enc{Härdle}{Haerdle}, W. K., Rychlik, T. (2010).
  \emph{Copula Theory and Its Applications}
  Springer, Lecture Notes in Statistics -- Proceedings.

  Schmid, F., Schmidt, R. (2007). Multivariate conditional versions of
  Spearman's rho and related measures of tail dependence.
  \emph{Journal of Multivariate Analysis} \bold{98}, 1123--1140.
  \doi{10.1016/j.jmva.2006.05.005}

  Schmidt, R., Stadtmueller, U. (2006). Non-parametric Estimation of
  Tail Dependence.
  \emph{Scandinavian Journal of Statistics} \bold{33}(2), 307--335.
  \doi{10.1111/j.1467-9469.2005.00483.x}
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
(lam.C <- c(lower = fitLambda(U, p = p)[2,1],
            upper = fitLambda(U, p = p, lower.tail = FALSE)[2,1])) # estimate lambdas
## => pretty good
U. <- pobs(U) # pseudo-observations (realistic)
(lam.C. <- c(lower = fitLambda(U., p = p)[2,1],
             upper = fitLambda(U., p = p, lower.tail = FALSE)[2,1])) # estimate lambdas
## => The pseudo-observations do have an effect...

## Higher-dimensional case
d <- 5
cop <- claytonCopula(2, dim = d)
set.seed(271)
U <- rCopula(n, copula = cop) # generate observations (unrealistic)
(lam.true <- lambda(cop)) # true tail-dependence coefficients lambda
(Lam.C <- list(lower = fitLambda(U, p = p),
               upper = fitLambda(U, p = p, lower.tail = FALSE))) # estimate Lambdas
## => Not too good
U. <- pobs(U) # pseudo-observations (realistic)
(Lam.C. <- list(lower = fitLambda(U., p = p),
                upper = fitLambda(U., p = p, lower.tail = FALSE))) # estimate Lambdas
## => Performance not too great here in either case
}
\keyword{nonparametric}
\keyword{multivariate}
