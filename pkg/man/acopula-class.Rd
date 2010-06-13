\name{acopula-class}
\Rdversion{1.1}
\title{Class "acopula" of Archimedean Copula Families}
\docType{class}
\alias{acopula-class}
%
\alias{initialize,acopula-method}
\alias{setTheta}
\alias{setTheta,acopula-method}
\alias{show,acopula-method}
%
\description{
  This class \code{"acopula"} of Archimedean Copula Families
  is mainly used for providing objects of known Archimedean families
  with all related functions.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("acopula", ...)}.
  For several well-known Archimedean copula families, the package \pkg{nacopula}
  already provides such family objects.
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{a string (class \code{"character"}) describing
      the copula family, e.g., "AMH" (or simply "A"), "Clayton" ("C"), "Frank" ("F"), "Gumbel" ("G"), or "Joe" ("J").}
    \item{\code{psi}, \code{psiInv}:}{The Archimedean generator, a \code{\link{function}}, and its inverse.}
    \item{\code{theta}:}{parameter value, a \code{\link{numeric}}, where
      \code{NA} means \dQuote{unspecified}.}
    \item{\code{paraInterval}:}{Either \code{\link{NULL}} or object of class
      \code{\linkS4class{"interval"}}, typically from a call such as
      \code{\link{interval}("[a,b)")}.}
    \item{\code{paraConstr}:}{A function of \code{theta} returning
      \code{TRUE} if and only if \code{theta} is a valid parameter value. Note
      that \code{paraConstr} is built automatically from the interval,
      whenever the \code{paraInterval} slot is valid.
      \code{\linkS4class{"interval"}}.}
	 \item{\code{nestConstr}:}{A \code{\link{function}}, which returns \code{TRUE} if and only if the two provided parameters \code{theta0} and \code{theta1} satisfy the sufficient nesting condition for this family.}
    \item{\code{V0}:}{A \code{\link{function}} which samples \code{n} random variates from the distribution \eqn{F} with Laplace-Stieltjes transform \eqn{\psi}{psi} and parameter \code{theta}.}
    \item{\code{V01}:}{A \code{\link{function}} which obtains a vector of realizations of \code{V0} and two parameters \code{theta0}, \code{theta1} which satisfy the sufficient nesting condition for this family, and which returns a vector of the same length as \code{V0} with random variates from the distribution function \eqn{F_{01}}{F01} with Laplace-Stieltjes transform \eqn{\psi_{01}}{psi01} and parameters \eqn{\vartheta_0=\,}{theta0=}\code{theta0}, \eqn{\vartheta_1=\,}{theta1=}\code{theta1}.}
    \item{\code{tau}, \code{tauInv}:}{Compute Kendall's tau of the bivariate Archimedean copula with generator \eqn{\psi}{psi} as a \code{\link{function}} of \code{theta}, respectively, \code{theta} as a function of Kendall's tau.}
    \item{\code{lambdaL} (\code{lambdaU}), \code{lambdaLInv} (\code{lambdaUInv}):}{Compute the lower (upper) tail-dependence coefficient of the bivariate Archimedean copula with generator \eqn{\psi}{psi} as a \code{\link{function}} of \code{theta}, respectively, \code{theta} as a function of the lower (upper) tail-dependence coefficient.}
  }
  For more details about Archimedean families, corresponding distributions and properties, see the references.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "acopula")}: is used to
      automatically construct the \code{paraConstr} function slot, when
      the \code{paraInterval} is provided (typically via
      \code{\link{interval}()}).}

    \item{setTheta}{\code{signature(x = "acopula")}: To set or change the
      parameter \eqn{\vartheta}{theta} (\code{theta}) of such an object, say
      \code{myCop}, use \code{setTheta(myCop, <value>)}.  This is used
      for construction of copula \dQuote{models} with a specified
      parameter.
    }
    \item{show}{\code{signature("acopula")}: compact overview of the copula.}
  }
}
\author{Martin Maechler}
\seealso{
  Specific provided copula family objects, e.g.,
  \code{\link{copAMH}}, \code{\link{copClayton}}, \code{\link{copFrank}},
  \code{\link{copGumbel}}, \code{\link{copJoe}}.

  A \emph{nested} Archimedean copula \emph{without}
  child copulas, see class \code{\linkS4class{"nacopula"}},
  is a proper Archimedean copula, and hence,
  \code{\link{onacopula}()} can be used to construct a specific
  parametrized Archimedean copula, see the example below.
}
\references{
	Hofert, M. (2010a),
	Efficiently sampling nested Archimedean copulas,
	\emph{Computational Statistics & Data Analysis}, in press.
	
	Hofert, M. (2010b),
	\emph{Sampling Nested Archimedean Copulas with Applications to CDO Pricing},
	Suedwestdeutscher Verlag fuer Hochschulschriften AG & Co. KG.
}
\examples{

## acopula class information
showClass("acopula")

## Information and structure of Clayton copulas
copClayton
str(copClayton)

## What are admissible parameters for Clayton copulas?
copClayton@paraInterval

## Can two Clayton copulas with parameters theta0 and theta1 be nested?
## Case 1: theta0 = 3, theta1 = 2 
theta0 <- 3; theta1 <- 2
copClayton@nestConstr(theta0,theta1) # FALSE
## Case 2: theta0 = 2, theta1 = 3 
theta0 <- 2; theta1 <- 3
copClayton@nestConstr(theta0,theta1) # TRUE

## Construct a bivariate Clayton copula with parameter theta
theta <- 2
C2 <- onacopula("Clayton", C(theta, 1:2))
C2@copula # is an "acopula" with specific parameter theta

## Plot the corresponding generator 
curve(C2@copula@psi(x, C2@copula@theta), xlim = c(0,5), ylim = c(0,1), col = 4)

## What is the corresponding Kendall's tau?
C2@copula@tau(theta)

## What are the corresponding tail-dependence coefficients?
C2@copula@lambdaL(theta)
C2@copula@lambdaU(theta)

## Generate n pairs of random variates from this copula
n <- 1000
U <- rnacopula(n, C2)

## Plot the generated pairs of random variates
splom2(U)

\keyword{classes}
