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
  This class \code{"acopula"} Archimedean Copula Families
  is mainly used for providing objects of known A.copula families
  with all related functions.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("acopula", ...)}.
  For half dozen known important copula families, the \pkg{nacopula}
  package already provides such family objects
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{a string (class \code{"character"}) describing
      the copula family.}
    \item{\code{psi}, \code{psiInv}:}{\eqn{psi}, the generator
      \code{\link{function}} and its inverse.}
    \item{\code{theta}:}{parameter value, a \code{\link{numeric}}, where
      \code{NA} means \dQuote{unspecified}.}
    \item{\code{paraInterval}:}{Either \code{\link{NULL}} or object of class
      \code{\linkS4class{"interval"}}, typically from a call such as
      \code{\link{interval}("[a,b)")}.}
    \item{\code{paraConstr}:}{A function of \code{theta} returning
      \code{TRUE} iff \code{theta} is a valid parameter value.  Note
      that \code{paraConstr} is built automatically from the interval,
      whenever the \code{paraInt} slot is a valid
      \code{\linkS4class{"interval"}}.}

%% TODO:
    \item{\code{V0}:}{Object of class \code{"function"} ~~ }
    \item{\code{tau}:}{Object of class \code{"function"} ~~ }
    \item{\code{tauInv}:}{Object of class \code{"function"} ~~ }
    \item{\code{lambdaL}:}{Object of class \code{"function"} ~~ }
    \item{\code{lambdaLInv}:}{Object of class \code{"function"} ~~ }
    \item{\code{lambdaU}:}{Object of class \code{"function"} ~~ }
    \item{\code{lambdaUInv}:}{Object of class \code{"function"} ~~ }
    \item{\code{nestConstr}:}{Object of class \code{"function"} ~~ }
    \item{\code{V01}:}{Object of class \code{"function"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "acopula")}: is used to
      automatically construct the \code{paraConstr} function slot, when
      the \code{paraInterval} is provided (typically via
      \code{\link{interval}()}).}

    \item{SetTheta}{\code{signature(x = "acopula")}: To set or change the
      parameter \eqn{\theta} (\code{theta}) of such an object, say
      \code{myCop}, use \code{setTheta(myCop, <value>)}.  This is used
      for construction of actual copula \dQuote{models} with a specified
      parameter.
    }
    \item{show}{\code{signature("acopula")}: compact overview of the copula.}
  }
}
% \author{Martin Maechler}
\seealso{
  Concrete provided copula family objects, such as
  \code{\link{copAMH}}, \code{\link{copClayton}}, etc.

  The special case of \emph{nested} Archimedean copula \emph{without}
  child copulas, see class \code{\linkS4class{"nacopula"}},
  is a (``proper'', unnested) Archimedean copula, and hence,
  \code{\link{onacopula}()} can be used to construct a concrete
  parametrized Archimedean copula, see the example below.
}
\examples{
showClass("acopula")
copClayton
str(copClayton)

C3 <- onacopula("Clayton", C(.25, 1:3))
C3@copula # is an "acopula" *with* specific parameter theta = 0.25

curve(C3@copula@psi(x, C3@copula@theta), 0, 3,
      n=1000, col=2, xaxs="i",yaxs="i")
}
\keyword{classes}
