\name{gnacopula}
\alias{gnacopula}
\alias{Rosenblatt}
\alias{gnacopulatrafo}
\title{Goodness-of-Fit Testing for (Nested) Archimedean Copulas}
\description{
  \code{gnacopula()} conducts a goodness-of-fit test for the given H0
  copula \code{cop} based on the data \code{u}.

  \code{Rosenblatt()} transforms vectors of random variates following
  the given (nested) Archimedean copula (with specified parameters) to
  \eqn{U[0,1]^d} vectors of random variates via the transformation of
  Rosenblatt (1952).

  \code{gnacopulatrafo()} transforms vectors of random variates following
  the given (nested) Archimedean copula (with specified parameters) to
  \eqn{U[0,1]^d} vectors of random variates via the transformation of
  Hering and Hofert (2011).
}
\usage{
gnacopula(u, cop, n.bootstrap = 0,
          estimation.method = eval(formals(enacopula)$method),
          include.K = TRUE, n.MC = 0, method = c("normal","log"),
          verbose = TRUE, \dots)
Rosenblatt(u, cop, n.MC = 0)
gnacopulatrafo(u, cop, include.K = TRUE, n.MC = 0)
}
\arguments{
  \item{u}{\eqn{n\times d}{n x d}-matrix of (pseudo-)observations (each
    value in \eqn{[0,1]}) from the copula to be tested.  Consider applying the
    function \code{\link{pobs}} first in order to obtain \code{u}.}
  \item{cop}{\eqn{H_0}{H0}-\code{"\linkS4class{outer_nacopula}"} with specified
    parameters to be tested for (currently only Archimedean copulas are provided).}
  \item{n.bootstrap}{non-negative integer giving the number of bootstrap
    replications. If \code{n.bootstrap=0}, then \code{\link{g01}} is
    applied to \code{\link{gnacopulatrafo}} with the chosen parameters.}
  \item{estimation.method}{estimation method, a \code{\link{character}}
    string used to determine the estimation method to be applied, see \code{\link{enacopula}}.}
  \item{include.K}{logical indicating whether the last component, involving
    \code{\link{K}}, is also used or not.}
  \item{n.MC}{parameter \code{n.MC} for \code{\link{K}}.}
  \item{method}{for \code{\link{g01}}, see there.}
  \item{\dots}{additional arguments to \code{\link{enacopula}}.}
  \item{verbose}{if TRUE, the progress of the bootstrap is displayed.}
}
\details{
  \describe{
    \item{\code{gnacopula}}{The bootstrap is explained in Genest et al. (2009), see also the source code.}

    \item{\code{Rosenblatt}}{Given a \eqn{d}-dimensional random vector
      \eqn{\bm{U}} following an Archimedean copula \eqn{C} with
      generator \eqn{\psi}, the conditional copula of \eqn{U_j=u_j} given
      \eqn{U_1=u_1,\dots,U_{j-1}=u_{j-1}}{U_1=u_1,...,U_{j-1}=u_{j-1}} is
      given by
      \deqn{C(u_j\,|\,u_1,\dots,u_{j-1})=\frac{\psi^{(j-1)}\Bigl(\sum_{k=1}^j\psi^{(-1)}(u_k)\Bigr)}{\psi^{(j-1)}\Bigl(\sum_{k=1}^{j-1}\psi^{(-1)}(u_k)\Bigr)}.}{%
	C(u_j | u_1,...,u_{j-1}) = (psi^{(j-1)}(sum(k=1..j )  psi^{(-1)}(u_k))) /
                                   (psi^{(j-1)}(sum(k=1..j-1) psi^{(-1)}(u_k))).}
      This formula is either evaluated with the exact derivatives or, if \code{n.MC} is positive, via Monte Carlo, see \code{\link{psiDabsMC}}.}

    \item{\code{gnacopulatrafo}}{Given a \eqn{d}-dimensional random vector \eqn{\bm{U}} following an
      Archimedean copula \eqn{C} with generator \eqn{\psi}, Hering and
      Hofert (2011) showed that \eqn{\bm{U}^\prime}{U'} is distributed
      as \eqn{\mathrm{U}[0,1]^d}{U[0,1]^d}, where
      \deqn{U_{j}^\prime=\left(\frac{\sum_{k=1}^{j}\psi^{-1}(U_{k})}{\sum_{k=1}^{j+1}\psi^{-1}(U_{k})}\right)^{j},\
	j\in\{1,\dots,d-1\},\ U_{d}^\prime=K(C(\bm{U})).}{%
	U'_j=((psi^{-1}(U_1)+...+psi^{-1}(U_j)) /
	      (psi^{-1}(U_1)+...+psi^{-1}(U_{j+1})))^j,  j in {1,...,d-1},
	U'_d=K(C(U)).}
      \code{gnacopulatrafo} applies this transformation row-wise to
      \code{u} and thus either a \eqn{n\times d}{n x d}-matrix or
      \eqn{n\times (d-1)}{n x (d-1)}-matrix, depending on whether the last
      component \eqn{U^\prime_d}{U'_d} which involves the possibly
      numerically challenging Kendall distribution function \eqn{K} is used
      (\code{include.K=TRUE}) or not (\code{include.K=FALSE})}
  }
}
\value{
  \code{gnacopula} returns an \R object of class \code{"htest"}.
  In the case \code{n.bootstrap = 0}, this object is simply the return value
  of \code{\link{g01}(\link{gnacopulatrafo}(...))}, otherwise the result also
  contains a list with the bootstrap results including the components
  \describe{
    \item{\code{AD.test}:}{the result of \code{\link{ad.test}} applied
      to \code{\link{g01trafo}(\link{gnacopulatrafo}(...))} (and thus
      the outcome of \code{\link{g01}(\link{gnacopulatrafo}(...))});}

    \item{\code{bootStats}:}{a list with component \code{theta.hat}
      containing the estimators of all bootstrap replications and
      component \code{AD.test} containing the results of
      \code{\link{ad.test}} for all bootstrap replications.}
  }

  \code{Rosenblatt} returns an \eqn{n\times d}{n x d}-matrix containing the transformed input \code{u}.

  \code{gnacopulatrafo} returns an
  \eqn{n\times d}{n x d}- (for \code{include.K=TRUE}) or
  \eqn{n\times (d-1)}{n x (d-1)}-matrix (for \code{include.K=FALSE})
  containing the transformed input \code{u}.
}
\author{Marius Hofert, Martin Maechler.}
\references{
  Genest, C., R\enc{é}{e}millard, B., and Beaudoin, D. (2009),
  Goodness-of-fit tests for copulas: A review and a power study
  \emph{Insurance: Mathematics and Economics}, \bold{44}, 199--213.

  Rosenblatt, M. (1952),
  Remarks on a Multivariate Transformation,
  \emph{The Annals of Mathematical Statistics}, \bold{23}, 3, 470--472.

  Hering, C. and Hofert, M. (2011),
  Goodness-of-fit tests for Archimedean copulas in large dimensions,
  submitted.

  Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2011b),
  Likelihood inference for Archimedean copulas,
  to be submitted.
}
\seealso{
  \code{\link{g01}} for the transformation(s) to the univariate setting.
}
\examples{
tau <- 0.5
(theta <- copGumbel@tauInv(tau)) # 2
d <- 5
(copG <- onacopulaL("Gumbel", list(theta,1:d)))

set.seed(1)
n <- 1000
X <- rnacopula(n, copG)
X <- qnorm(X) # x now follows a meta-Gumbel model with N(0,1) marginals
U <- pobs(X) # build pseudo-observations

## check if the data comes from a meta-Gumbel model
## (choose larger n.bootstrap in a realistic setup)
res.H0.G <- gnacopula(U, cop=copG, n.bootstrap=10,
                      estimation.method="mle")
res.H0.G$p.value # non-rejection according to 5\% level

## plot of the transformed data (Rosenblatt (1952)) to check uniformity
U. <- Rosenblatt(U, cop=copG) # exact
pairs(U., cex=0.2) # looks good

## plot of the transformed data (Rosenblatt (1952)) to check uniformity
U. <- Rosenblatt(U, cop=copG, n.MC=10000) # via Monte Carlo
pairs(U., cex=0.2) # looks good

## plot of the transformed data (Hering and Hofert (2011)) to check uniformity
U. <- gnacopulatrafo(U, cop=copG)
pairs(U., cex=0.2) # looks good

## what about a meta-Clayton model?
## (choose larger n.bootstrap in a realistic setup)
## note: the parameter of the Clayton copula is only a dummy,
##       it will be estimated anyway
copC <- onacopulaL("Clayton", list(1,1:d))
res.H0.C <- gnacopula(U, cop=copC, n.bootstrap=10,
	              estimation.method="mle")
res.H0.C$p.value # rejection according to 5\% level

## plot of the transformed data (Rosenblatt (1952)) to see the
## deviations from uniformity
U. <- Rosenblatt(U, cop=copC, n.MC=0)
pairs(U., cex=0.2) # clearly visible

## plot of the transformed data (Hering and Hofert (2011)) to see the
## deviations from uniformity
U. <- gnacopulatrafo(U, cop=copC, include.K=TRUE, n.MC=0)
pairs(U., cex=0.2) # clearly visible

## and for the map to a supposedly U[0,1] distribution?
Z <- g01trafo(U.)
plot(1:length(Z), Z)
#-> harder to see, but still does not look perfectly uniform (as expected)
}
\keyword{htest}
\keyword{distribution}
\keyword{multivariate}

