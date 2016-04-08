\name{fitCopula}
\alias{fitCopula}
\alias{fitCopula-methods}
\alias{fitCopula,copula-method}
\alias{fitCopula,rotCopula-method}
\alias{loglikCopula}
\title{Copula Parameter Estimation}
\description{
  Parameter estimation (fitting) of a copula model to multivariate data.
}
\usage{
loglikCopula(param, u, copula)
fitCopula(copula, data, method = c("mpl", "ml", "itau", "irho", "itau.mpl"),
          posDef = is(copula, "ellipCopula"),
          start = NULL, lower = NULL, upper = NULL,
          optim.method = "BFGS", optim.control = list(maxit=1000),
          estimate.variance = NA, hideWarnings = FALSE, \dots)
}
\arguments{
  \item{param}{A vector of parameter values.}
  \item{u}{\eqn{n\times d}{n x d}-matrix of (pseudo-)observations in
    \eqn{[0,1]^d} for computing the copula log-likelihood, where \eqn{n} denotes
    the sample size and \eqn{d} the dimension. Consider applying the function
    \code{\link{pobs}()} first in order to obtain such data.}
  \item{data}{As \eqn{u}, an \eqn{n\times d}{n x d}-matrix of data. For
    \code{method} being \code{"mpl"}, \code{"ml"} or \code{"itau.mpl"},
    this has to be data in \eqn{[0,1]^d}. For \code{method} being
    \code{"itau"} or \code{"irho"}, it can either be data in \eqn{[0,1]^d}
    or in the whole \eqn{d}-dimensional space.}
  \item{copula}{A \code{"\linkS4class{copula}"} object.}
  \item{method}{A \code{\link{character}} string specifying the copula
    parameter estimator used. This can be one of:
    \describe{
      \item{"mpl"}{Maximum pseudo-likelihood estimator (based on the given data
	in \eqn{[0,1]^d}; typically obtained via \code{\link{pobs}()}).}
      \item{"ml"}{As \code{"mpl"} just with a different variance
	estimator. For this to be correct (thus giving the true MLE),
	\code{data} must contain observations from the true
	underlying copula to be estimated (assumed to be known in this case).}
      \item{"itau"}{Inversion of Kendall's tau estimator. \code{data}
	can be either in \eqn{[0,1]^d} (true or pseudo-observations of
	the underlying copula to be estimated) or in the \eqn{d}-dimensional space.}
      \item{"irho"}{As \code{"itau"} just with Spearman's rho instead of
	Kendall's tau.}
      \item{"itau.mpl"}{This is the estimator of \eqn{t} copula
	parameters suggested by Demarta and McNeil (2005) based on the
	given \code{data} in \eqn{[0,1]^d} (true or pseudo-observations of
	the underlying copula to be estimated).}
    }
  }
  \item{posDef}{A \code{\link{logical}} indicating whether a proper
    correlation matrix is computed.}
  \item{start}{A \code{\link{vector}} of starting values for the
    parameter optimization via \code{\link{optim}()}.}
  \item{lower, upper}{Lower or Upper parameter bounds for the
    optimization methods \code{"Brent"} or \code{"L-BFGS-B"}.}
  \item{optim.control}{A \code{\link{list}} of control parameters
    passed to \code{\link{optim}(*, control=optim.control)}.}
  \item{optim.method}{The optimization method; see
    \code{\link{optim}()}'s \code{method}.}
  \item{estimate.variance}{A \code{\link{logical}} indicating whether
    the estimator's asymptotic variance is computed (if available for
    the given \code{copula}; the default \code{NA} computes it for the
    \code{method}s \code{"itau"} and \code{"irho"}, cannot (yet) compute
    it for \code{"itau.mpl"} and only computes it for \code{"mpl"} or \code{"ml"}
    if the optimization converged).}
  \item{hideWarnings}{A \code{\link{logical}}, which, if
    \code{\link{TRUE}}, suppresses warnings from the involved
    likelihood maximization (typically when the likelihood is
    evaluated at invalid parameter values).}
  \item{\dots}{Additional arguments passed to the underlying auxiliary functions.}
}
\value{
  \code{loglikCopula()} returns the copula log-likelihood evaluated at the
  parameter (vector) \code{param} given the data \code{u}.

  The return value of \code{fitCopula()} is an object of class
  \code{"\linkS4class{fitCopula}"} (see there), containing (among
  others!) the slots
  \describe{
    \item{estimate}{The parameter estimates.}
    \item{var.est}{The large-sample (i.e., asymptotic) variance estimate of the parameter
      estimator; (filled with) \code{NA} if \code{estimate.variance=FALSE}.
    }
    \item{copula}{The fitted copula object.}
  }

  The \code{\link{summary}()} method for \code{"fitCopula"} objects
  returns an S3 \dQuote{class} \code{"summary.fitCopula"}, which is simply a list
  with components \code{method}, \code{loglik} and \code{convergence},
  all three from the corresponding slots of the
  \code{"\linkS4class{fitCopula}"} objects, and \code{coefficients}
  (a matrix of estimated coefficients, standard errors, t values and p-values).
}
\details{
  The only difference between \code{"mpl"} and \code{"ml"} is in the
  variance-covariance estimate, \emph{not} in the parameter
  (\eqn{\theta}) estimates.

  If method \code{"mpl"} in \code{fitCopula()} is used and if
  \code{start} is not assigned a value, estimates obtained
  from method \code{"itau"} are used as initial values in
  the optimization.

  In the multiparameter elliptical case and when the estimation is based
  on Kendall's tau or Spearman's rho, the estimated correlation matrix
  may not always be positive-definite. In that case,
  \code{\link[Matrix]{nearPD}(*, corr=TRUE)} (from \CRANpkg{Matrix}) is
  applied to get a proper correlation matrix.

  If methods \code{"itau"} or \code{"irho"} are used in
  \code{fitCopula()}, an estimate of the asymptotic variance (if
  available for the copula under consideration) will be correctly
  computed only if the argument \code{data} consists of pseudo-observations
  (see \code{\link{pobs}()}).

  Consider the \eqn{t} copula with \code{df.fixed=FALSE} (see \code{\link{ellipCopula}()}).
  In this case, the methods \code{"itau"} and \code{"irho"} cannot be used in
  \code{fitCopula()} as they cannot estimate the degrees of freedom
  parameter \code{df}. For the methods \code{"mpl"} and \code{"itau.mpl"}
  the asymptotic variance cannot be estimated (yet). And for the methods
  \code{"ml"} and \code{"mpl"}, when \code{start} is not specified, the
  starting value for \code{df} is set to \code{copula@df}, typically 4.

  To implement the \emph{Inference Functions for Margins (IFM} method
  (see, e.g., Joe 2005), set \code{method="ml"} and note that
  \code{data} need to be pseudo-observations obtained from \emph{fitted}
  parametric marginal distribution functions. The returned large-sample
  variance will then underestimate the true variance (as the procedure
  cannot take into account the (unknown) estimation error for the margins).

  The fitting procedures based on \code{\link{optim}()} generate
  warnings because invalid parameter values are tried during the
  optimization process. When the number of parameters is one and the
  parameter space is bounded, using \code{optim.method="Brent"} is
  likely to give less warnings. Furthermore, from experience,
  \code{optim.method="Nelder-Mead"} is sometimes a more robust
  alternative to \code{optim.method="BFGS"}.

  There are methods for \code{\link{vcov}()}, \code{\link{coef}()},
  \code{\link{logLik}()}, and \code{\link{nobs}()}.
}
%\author{Marius Hofert, Ivan Kojadinovic, Martin Maechler and Jun Yan}
\references{
  Genest, C. (1987). Frank's family of bivariate distributions.
  \emph{Biometrika} \bold{74}, 549--555.

  Genest, C. and Rivest, L.-P. (1993). Statistical inference procedures
  for bivariate Archimedean copulas.
  \emph{Journal of the American Statistical Association} \bold{88}, 1034--1043.

  Rousseeuw, P. and Molenberghs, G. (1993). Transformation of nonpositive
  semidefinite correlation matrices.
  \emph{Communications in Statistics: Theory and Methods} \bold{22}, 965--984.

  Genest, C., Ghoudi, K., and Rivest, L.-P. (1995). A semiparametric
  estimation procedure of dependence parameters in multivariate
  families of distributions. \emph{Biometrika} \bold{82}, 543--552.

  Joe, H. (2005). Asymptotic efficiency of the two-stage estimation
  method for copula-based models.
  \emph{Journal of Multivariate Analysis} \bold{94}, 401--419.

  Demarta, S. and McNeil, A. J. (2005). The t copula and related copulas.
  \emph{International Statistical Review} \bold{73}, 111--129.

  Genest, C. and Favre, A.-C. (2007). Everything you always wanted to know about
  copula modeling but were afraid to ask.
  \emph{Journal of Hydrologic Engineering} \bold{12}, 347--368.

  Kojadinovic, I. and Yan, J. (2010). Comparison of three semiparametric methods
  for estimating dependence parameters in copula models.
  \emph{Insurance: Mathematics and Economics} \bold{47}, 52--63.
}
\seealso{
  \code{\link{Copula}},
  \code{\link{mvdc}} for fitting multivariate distributions
  \emph{including} the margins aka \dQuote{meta copula}s;
  \code{\link{gofCopula}}.

  For maximum likelihood of (nested) Archimedean copulas, see
  \code{\link{emle}}, etc.
}
\examples{
(Xtras <- copula:::doExtras()) # determine whether examples will be extra (long)
n <- if(Xtras) 200 else 64 # sample size

## A Gumbel copula
set.seed(7) # for reproducibility
gumbel.cop <- gumbelCopula(3, dim=2)
x <- rCopula(n, gumbel.cop) # "true" observations (simulated)
u <- pobs(x)                # pseudo-observations
## Inverting Kendall's tau
fit.tau <- fitCopula(gumbel.cop, u, method="itau")
fit.tau
coef(fit.tau)# named vector
## Inverting Spearman's rho
fit.rho <- fitCopula(gumbel.cop, u, method="irho")
fit.rho
## Maximum pseudo-likelihood
fit.mpl <- fitCopula(gumbel.cop, u, method="mpl")
fit.mpl
## Maximum likelihood
fit.ml <- fitCopula(gumbel.cop, x, method="ml")
fit.ml # print()ing works via summary() ...
## ... and of that, what's the log likelihood (in two different ways):
(ll. <- logLik(fit.ml))
stopifnot(all.equal(as.numeric(ll.),
    loglikCopula(coef(fit.ml), u=x, copula=gumbel.cop)))

## A Gauss/normal copula

## With multiple/*un*constrained parameters
set.seed(6) # for reproducibility
normal.cop <- normalCopula(c(0.6, 0.36, 0.6), dim=3, dispstr="un")
x <- rCopula(n, normal.cop) # "true" observations (simulated)
u <- pobs(x)                # pseudo-observations
## Inverting Kendall's tau
fit.tau <- fitCopula(normal.cop, u, method="itau")
fit.tau
## Inverting Spearman's rho
fit.rho <- fitCopula(normal.cop, u, method="irho")
fit.rho
## Maximum pseudo-likelihood
fit.mpl <- fitCopula(normal.cop, u, method="mpl")
fit.mpl
coef(fit.mpl) # named vector
str(sf.mpl <- summary(fit.mpl))
coef(sf.mpl) # the matrix, with SE, t-value, ...
## Maximum likelihood
fit.ml <- fitCopula(normal.cop, x, method="ml")
fit.ml

## With dispstr="toep"
normal.cop.toep <- normalCopula(c(0, 0), dim=3, dispstr="toep")
## Inverting Kendall's tau
fit.tau <- fitCopula(normal.cop.toep, u, method="itau")
fit.tau
## Inverting Spearman's rho
fit.rho <- fitCopula(normal.cop.toep, u, method="irho")
fit.rho
## Maximum pseudo-likelihood
fit.mpl <- fitCopula(normal.cop.toep, u, method="mpl")
fit.mpl
## Maximum likelihood
fit.ml <- fitCopula(normal.cop.toep, x, method="ml")
fit.ml

## With dispstr="ar1"
normal.cop.ar1 <- normalCopula(c(0), dim=3, dispstr="ar1")
## Inverting Kendall's tau
fit.tau <- fitCopula(normal.cop.ar1, u, method="itau")
fit.tau
## Inverting Spearman's rho
fit.rho <- fitCopula(normal.cop.ar1, u, method="irho")
fit.rho
## Maximum pseudo-likelihood
fit.mpl <- fitCopula(normal.cop.ar1, u, method="mpl")
fit.mpl
## Maximum likelihood
fit.ml <- fitCopula(normal.cop.ar1, x, method="ml")
fit.ml

## A t copula with variable df (df.fixed=FALSE)
(tCop <- tCopula(c(0.2,0.4,0.6), dim=3, dispstr="un", df=5))
set.seed(101)
x <- rCopula(n, tCop) # "true" observations (simulated)
u <- pobs(x)          # pseudo-observations
## Maximum likelihood (start = (rho[1:3], df))
(tc.ml  <- fitCopula(tCop, x, method="ml", start=c(0,0,0, 10)))
(tc.ml. <- fitCopula(tCop, x, method="ml")) # without 'start'
## Maximum pseudo-likelihood (the asymptotic variance cannot be estimated)
(tc.mpl <- fitCopula(tCop, u, method="mpl", estimate.variance=FALSE,
                     start= c(0,0,0,10)))

## Extra checks
if(Xtras) { ## Typically not run on R CMD check

## Without 'start'
(tc.mp. <- fitCopula(tCop, u, method="mpl", estimate.variance=FALSE))
all.eqCop <- function(x,y, ...) {
   x@fitting.stats$counts <- y@fitting.stats$counts <- NULL
   all.equal(x,y, ...) }
stopifnot(all.eqCop(tc.ml , tc.ml., tolerance= .005),
          all.eqCop(tc.mpl, tc.mp., tolerance= .005))

## The same t copula but with df.fixed=TRUE (=> use the same data!)
tC.f <- tCopula(dim=3, dispstr="un", df=5, df.fixed=TRUE)
## Maximum likelihood (start = rho[1:3])
(tcF.ml  <- fitCopula(tC.f, x, method="ml", start=c(0,0,0)))
(tcF.ml. <- fitCopula(tC.f, x, method="ml"))# without 'start'
stopifnot(all.eqCop(tcF.ml,tcF.ml., tolerance= 4e-4))
vcov(tcF.ml) # the (estimated, asymptotic) var-cov matrix

## Maximum pseudo-likelihood (the asymptotic variance cannot be estimated)
(tcF.mpl <- fitCopula(tC.f, u, method="mpl", estimate.variance=FALSE,
                      start=c(0,0,0)))
(tcF.mp. <- fitCopula(tC.f, u, method="mpl", estimate.variance=FALSE))
stopifnot(all.eqCop(tcF.mpl,tcF.mp., tolerance= 1e-5))

} # end Xtras
}
\keyword{models}
\keyword{multivariate}
