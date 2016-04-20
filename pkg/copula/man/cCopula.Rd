 \name{cCopula}
\alias{cCopula}
\alias{cacopula}
\title{Conditional Copula Function}
\description{
  Compute the conditional distribution aka conditional copula function,
  \eqn{C(u_d\,|\,u_1,\dots, u_{d-1})}{C(u[d] | u[1],..., u[d-1])}
  of \eqn{u_d}{u[d]} given \eqn{u_1,\dots,u_{d-1}}{u[1],..., u[d-1]}.
}
\usage{
cCopula(u, copula, j.ind=ncol(u), n.MC=0, inverse=FALSE, log=FALSE)

## Deprecated, for Archimedean copulas only:
cacopula(u, cop, n.MC=0, log=FALSE)
}
\arguments{
  \item{u}{\eqn{n\times d}{n x d}-matrix; the conditioning is done on the values
    in the first \eqn{d-1} columns.}
  \item{copula}{any \code{"\linkS4class{copula}"} with specified
    parameters; currently only Archimedean and elliptical copulas
    provided.}
  \item{j.ind}{indices \eqn{j\ge 2}{j>=2} for which
    \eqn{C(u_j\,|\,u_1,\dots,u_{j-1})}{C(u[j]|u[1],...,u[j-1])} is computed.}
  \item{cop}{\code{"\linkS4class{outer_nacopula}"} with specified
    parameters (only Archimedean copulas are currently provided).}
  \item{n.MC}{Monte Carlo sample size; for Archimedean copulas only,
    used if positive.}
  \item{inverse}{logical indicating whether the inverse
    conditional copula is returned.}
  \item{log}{if TRUE the logarithm of the conditional copula function is
    returned.}
}
\value{
  \code{\link{numeric}} vector of length \eqn{n} containing the
  conditional copula function of \eqn{u_d}{u[d]} given
  \eqn{u_1,\dots,u_{d-1}}{u[1],...,u[d-1]}.
}
\details{
  \code{cCopula()} can be seen as a special case of the Rosenblatt
  transform which by default computes all conditional distributions.
  Conseqently, \code{cCopula()} simply calls \code{\link{rtrafo}()}.

  Note that \emph{sampling} aka random number generation from the
  conditional distribution can be achieved by using \code{inverse=TRUE}.
}
%\author{Marius Hofert, Yongsheng Wang, and Martin Maechler}
\note{
  For some (but not all) families, this function also makes sense on the
  boundaries (if the corresponding limits can be computed).
}
\seealso{
  \code{\link{rtrafo}}; \code{\link{acopula-families}}.
}
\examples{
## Note: This illustrates the numerically more involved (and thus much slower)
##       Gumbel case.

## Define the copula
tau <- 0.5
theta <- iTau(gumbelCopula(), tau = tau)
d <- 5
gc <- gumbelCopula(theta, dim = d)

## Sample via conditional distribution method
n <- 200
set.seed(271)
U. <- matrix(runif(n*d), ncol = d) # U(0,1)^d
if(FALSE)
    Rprof(profiling <- tempfile(), line.profiling=TRUE) # enable profiling
U <- cCopula(U., copula = gc, j.ind = NULL, inverse = TRUE) # slow for ACs except Clayton
if(FALSE) {
    Rprof(NULL) # disable profiling
    (profile <- summaryRprof(profiling, lines="both")) # get a summary
}
splom2(U) # scatter-plot matrix copula sample

## Rosenblatt-transform back to U(0,1)^d
U. <- cCopula(U, copula = gc, j.ind = NULL)
splom2(U.) # U(0,1)^d again
}
\keyword{distribution}
