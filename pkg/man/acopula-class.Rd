\name{acopula-class}
\Rdversion{1.1}
\title{Class "acopula" of Archimedean Copula Families}
\docType{class}
\alias{acopula-class}
%
\alias{initialize,acopula-method}
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
    \item{\code{name}:}{A string (class \code{"character"}) describing
      the copula family, e.g., "AMH" (or simply "A"), "Clayton" ("C"),
      "Frank" ("F"), "Gumbel" ("G"), or "Joe" ("J").}
    \item{\code{psi}, \code{psiInv}:}{The Archimedean generator \eqn{\psi}{psi}
 	   and its inverse (\code{\link{function}}).}
    \item{\code{psiDAbs}:}{A \code{\link{function}} which computes the absolute
      value of the derivative of the generator \eqn{\psi}{psi} for the given
 		parameter \code{theta} and of the given \code{degree} (defaults to 1).
 		Note that there is not loss in information since the derivatives alternate
 		in sign so that the generator derivative itself can be accessed with
 		(-1)^\code{degree}\code{psiDAbs}. If 
 		the boolean parameter \code{MC} is \code{FALSE} (the default) a series
 		representation for the derivative is used; otherwise, the derivative is
 		evaluated by Monte Carlo simulation. The parameter \code{N} denotes
 		the truncation point for the series representation if \code{MC} is
      \code{FALSE} (note that a good choice of \code{N} depends on \code{theta})
 		and the number of variates to be drawn if \code{MC} is \code{TRUE}. Both
 		\code{MC} and \code{N} do not exist for Clayton's Archimedean family since
 		all derivatives can be given explicitly
 		for this Archimedean family. For all families, the optional parameter
 		\code{log} (defaults to \code{FALSE}) indicates whether the logarithm of
 		the absolute value of the derivative is returned. Note that for
 		\code{degree} equal to 1, all	computations are exact.}
    \item{\code{theta}:}{Parameter value, a \code{\link{numeric}}, where
      \code{NA} means \dQuote{unspecified}.}
    \item{\code{paraInterval}:}{Either \code{\link{NULL}} or object of class
      \code{"\linkS4class{interval}"}, typically from a call such as
      \code{\link{interval}("[a,b)")}.}
    \item{\code{paraConstr}:}{A function of \code{theta} returning
      \code{TRUE} if and only if \code{theta} is a valid parameter value.  Note
      that \code{paraConstr} is built automatically from the interval,
      whenever the \code{paraInterval} slot is valid.
      \code{"\linkS4class{interval}"}.}
    \item{\code{nestConstr}:}{A \code{\link{function}}, which returns
      \code{TRUE} if and only if the two provided parameters \code{theta0} and
      \code{theta1} satisfy the sufficient nesting condition for this family.}
    \item{\code{V0}:}{A \code{\link{function}} which samples \code{n} random
      variates from the distribution \eqn{F} with Laplace-Stieltjes transform
      \eqn{\psi}{psi} and parameter \code{theta}.}
	 \item{\code{dV0}:}{A \code{\link{function}} which computes either the
	   probability mass function or the probability density function (depending
	 	on the Archimedean family) of the distribution function whose
	 	Laplace-Stieltjes transform equals the generator \eqn{\psi}{psi} at the
	 	argument \code{x} (possibly a vector) for the given parameter
	 	\code{theta}. An optional argument \code{log} indicates whether
	   or not the logarithm of the mass or density is computed (defaults to
	   \code{FALSE}).}
    \item{\code{V01}:}{A \code{\link{function}} which gets a vector of
      realizations of \code{V0} and two parameters \code{theta0}, \code{theta1}
      which satisfy the sufficient nesting condition for this family, and which
      returns a vector of the same length as \code{V0} with random variates from
      the distribution function \eqn{F_{01}}{F01} with Laplace-Stieltjes
      transform \eqn{\psi_{01}}{psi01} and parameters
      \eqn{\theta_0=\,}{theta0=}\code{theta0},
      \eqn{\theta_1=\,}{theta1=}\code{theta1}.}
	 \item{\code{dV01}:}{Similar to \code{dV0} with the difference being that
	 	this \code{\link{function}} computes the probability mass or density
	 	function for the Laplace-Stieltjes transform
 		\deqn{\psi_{01}(t;V_0)=\exp(-V_0\psi_0^{-1}(\psi_1(t)))}{psi01(t;V0)=exp(-V0*psi0^{-1}(psi1(t)))} corresponding to the distribution function \eqn{F_{01}}{F01}.
	
	   Arguments are the evaluation point(s) \code{x}, the value(s) \code{V0},
	   and the parameters \code{theta0} and \code{theta1}. As for \code{dV0}, the
	   optional argument \code{log} can be specified (defaults to \code{FALSE}).
	   Note that if \code{x} is a vector, \code{V0} must either have
	   length one (in which case \code{V0} is the same for every component of
	   \code{x}) or \code{V0} must be of the same length as \code{x} (in which
	   case the components of \code{V0} correspond to the ones of \code{x}).}
    \item{\code{cCdf}:}{A \code{\link{function}} with arguments \code{v},
 		\code{u}, and \code{theta} which computes the bivariate conditional copula
		function at \code{v} given \code{u}, i.e., \eqn{C(v\vert u)}{C(v|u)}.}
	 \item{\code{diag}:}{The copula diagonal.}
	 \item{\code{dDiag}:}{The density of the copula diagonal (note that the
	 	diagonal of a copula is a cumulative distribution function).}
 	 \item{\code{dAc}:}{A \code{\link{function}} which computes the
      logarithm of the density \eqn{c} for the Archimedean family under
      consideration, given by
		\deqn{\log c(u_1,\dots,u_d)=\log\frac{\psi^{(d)}(\sum_{j=1}^d\psi^{-1}(u_j))}{\prod_{j=1}^d\psi^{(d)}(\psi^{-1}(u_j))}.}{log(c(u))=log(psi^{(d)}(sum_{j=1}^d psi^{-1}(u_j))/prod_{j=1}^d psi^prime(psi^{-1}(u_j))).} Parameters of
		\code{dAc} are the \code{\link{matrix}} \code{u} (the number of
		columns of \code{u} determines \eqn{d}), the parameter \code{theta}, and
		(for the Archimedean family not being Clayton) the arguments
		\code{MC} and \code{N} (with the same meaning as for \code{psiDAbs}).}
 	\item{\code{K}:}{Computes the Kendall distribution \code{\link{function}} at
      \code{t} given the parameter \code{theta} and the dimension \code{d}. For
      an Archimedean copula with generator \eqn{\psi}, this function is given by
  \deqn{K(t)=\sum_{k=0}^{d-1}\frac{(-1)^k\psi^{(k)}(\psi^{-1}(t))}{k!}(\psi^{-1}(t))^k.}{K(t)=sum_{k=0}^{d-1}(-1)^k*psi^{(k)}((psi^{-1}(t))/k!)(psi^{-1}(t))^k.}
 		For all except Clayton's family, the arguments \code{MC} and
 		\code{N} can be specified (with the same meaning as for \code{psiDAbs}).}
    \item{\code{tau}, \code{tauInv}:}{Compute Kendall's tau of the bivariate
      Archimedean copula with generator \eqn{\psi}{psi} as a
      \code{\link{function}} of \code{theta}, respectively, \code{theta} as a
      function of Kendall's tau.}
    \item{\code{lambdaL}, \code{lambdaU}, \code{lambdaLInv},
 		\code{lambdaUInv}:}{Compute the lower (upper) tail-dependence coefficient
 		of the bivariate Archimedean copula with generator \eqn{\psi}{psi} as a
      \code{\link{function}} of \code{theta}, respectively, \code{theta} as a
      function of the lower (upper) tail-dependence coefficient.}
  }
  For more details about Archimedean families, corresponding distributions and
  properties, see the references.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "acopula")}: is used to
      automatically construct the \code{paraConstr} function slot, when
      the \code{paraInterval} is provided (typically via
      \code{\link{interval}()}).}

    \item{show}{\code{signature("acopula")}: compact overview of the copula.}
  }
}
\author{Martin Maechler, Marius Hofert}
\seealso{
  Specific provided copula family objects, e.g.,
  \code{\link{copAMH}}, \code{\link{copClayton}}, \code{\link{copFrank}},
  \code{\link{copGumbel}}, \code{\link{copJoe}}.
  \cr To access these use \code{\link{getAcop}}.

  A \emph{nested} Archimedean copula \emph{without}
  child copulas, see class \code{"\linkS4class{nacopula}"},
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
copClayton@nestConstr(theta0 = 3, theta1 = 2) # FALSE
## Case 2: theta0 = 2, theta1 = 3
copClayton@nestConstr(theta0 = 2, theta1 = 3) # TRUE

## For more examples, see  help("acopula-families")
}
\keyword{classes}
