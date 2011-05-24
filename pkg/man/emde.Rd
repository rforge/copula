\name{emde}
\alias{emde}
\title{Minimum Distance Estimators for (Nested) Archimedean Copulas}
\description{
  Compute minimum distance estimators for (nested) Archimedean copulas.
}
\usage{
emde(u, cop,
     method=c("mde.normal.CvM", "mde.normal.KS",
              "mde.log.CvM",    "mde.log.KS"),
     interval=paraOptInterval(u, cop@copula@name),
     include.K = ncol(u) <= 5, \dots)
}
\arguments{
  \item{u}{\eqn{n\times d}{n x d}-matrix of (pseudo-)observations (each value in \eqn{[0,1]}) from the copula, where \eqn{n} denotes the sample size and \eqn{d} the dimension.}
  \item{cop}{\code{\linkS4class{outer_nacopula}} to be estimated (currently only Archimedean copulas are provided).}
  \item{method}{a \code{\link{character}} string specifying the distance
    method, which has to be one (or a unique abbreviation) of
    \describe{
      \item{\code{"mde.normal.CvM"}}{map to a chi-square distribution and then
            use a Cramer-von Mises distance;}
      \item{\code{"mde.normal.KS"}}{map to a chi-square distribution and then
            use Kolmogorov-Smirnov distance;}
      \item{\code{"mde.log.CvM"}}{map to an Erlang distribution and then use a
            Cramer-von Mises distance;}
      \item{\code{"mde.log.KS"}}{map to an Erlang distribution and then use a
            Kolmogorov-Smirnov distance.}
    }
    The four methods are described in Hofert et al. (2011a).  See also the
    \sQuote{Details} section.}
  \item{interval}{bivariate vector denoting the interval where
    optimization takes place.  The default is computed as described in Hofert et
    al. (2011a).}
  \item{include.K}{logical indicating whether the last component, the Kendall
    distribution function K, is also used or not (see the \sQuote{Details}
    section for more information).}
  \item{\dots}{additional arguments passed to \code{\link{optimize}}.}
}
\details{
   Given a \eqn{d}-dimensional random vector \eqn{\bm{U}} following an
   Archimedean copula \eqn{C} with generator \eqn{\psi}, Hering and
   Hofert (2011) showed that \eqn{\bm{U}^\prime}{U^prime} is distributed
   as \eqn{\mathrm{U}[0,1]^d}{U[0,1]^d}, where
   \deqn{U_{j}^\prime=\left(\frac{\sum_{k=1}^{j}\psi^{-1}(U_{k})}{\sum_{k=1}^{j+1}\psi^{-1}(U_{k})}\right)^{j},\
     j\in\{1,\dots,d-1\},\ U_{d}^\prime=K(C(\bm{U})).}{%
     U_j^prime=((psi^{-1}(U_1)+...+psi^{-1}(U_j)) / (psi^{-1}(U_1)+...+psi^{-1}(U_{j+1})))^j, j in {1,...,d-1},
     U_d^prime=K(C(U)).}
   This transformation is first used to map the \eqn{n\times d}{n x d}-matrix
   of given realizations to a \eqn{n\times d}{n x d}-matrix or
   \eqn{n\times (d-1)}{n x (d-1)}-matrix, depending on whether the last
   component \eqn{U^\prime_d}{U^prime_d} which involves the possibly
   numerically challenging Kendall distribution function \eqn{K} is used
   (\code{include.K=TRUE}) or not (\code{include.K=FALSE}).  Then, using
   either the sum of the squares of the quantile function of the
   standard normal distribution (for \code{method="mde.normal.CvM"} and
   \code{method="mde.normal.KS"}) or the sum of the negative logarithms
   (for \code{method="mde.log.CvM"} and \code{method="mde.log.KS"}), the
   transformed quantities are mapped to a chi-square or an Erlang
   distribution, respectively.  Afterwards, a Cramér-von Mises (for
   \code{method="mde.normal.CvM"} and \code{method="mde.log.CvM"}) or
   Kolmogorov-Smirnov (for \code{method="mde.normal.KS"} and
   \code{method="mde.log.KS"}) distance is applied.  This is repeated in
   an optimization until the copula parameter is found such that this
   distance is minimized.

 Note that the same transformations as described above are applied for goodness-of-fit testing.
}
\value{
  \code{\link{list}} as returned by \code{\link{optimize}}, including the
  minimum distance estimator.
}
\author{Marius Hofert}
\references{
   Hofert, M., Mächler, M., and McNeil, A. J. (2011a),
   \emph{Estimators for Archimedean copulas in high dimensions: A comparison},
   to be submitted.

   Hering ...  Hofert %% FIXME
}
\seealso{
  \code{\link{enacopula}} (wrapper for different estimators),
  \code{\link{gnacopula}} (transformation \eqn{\bm{U}}{U} to
  \eqn{\bm{U}^\prime}{U^prime}),
  \code{\link{g01}} (transformation to a chi-square or Erlang
  distribution), and
  \code{\link{K}} (Kendall distribution function).
}
\examples{
n <- 200
d <- 20
tau <- 0.25
(theta <- copGumbel@tauInv(tau)) # 4/3
(cop <- onacopulaL("Gumbel", list(theta,1:d)))

set.seed(1)
U <- rnacopula(n, cop)

(meths <- eval(formals(emde)$method)) # "mde.normal.CvM", "mde.normal.KS", ...
fun <- function(meth, u, cop, theta){
	run.time <- system.time(val <- emde(u, cop=cop, method=meth)$minimum)
	list(value=val, error=val-theta, utime.ms=1000*run.time[[1]])
}
(res <- sapply(meths, fun, u=U, cop=cop, theta=theta))
}
\keyword{models}
