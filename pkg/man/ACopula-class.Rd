\name{acopula-class}
\Rdversion{1.1}
\title{Class "acopula" of Archimedean Copula Families}
\docType{class}
\alias{acopula-class}
%
\alias{initialize,acopula-method}
\alias{setTheta}
\alias{setTheta,acopula-method}
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
    \item{\code{lTDC}:}{Object of class \code{"function"} ~~ }
    \item{\code{lTDCInv}:}{Object of class \code{"function"} ~~ }
    \item{\code{uTDC}:}{Object of class \code{"function"} ~~ }
    \item{\code{uTDCInv}:}{Object of class \code{"function"} ~~ }
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
  }
}
% \references{ ~put references to the literature/web site here ~ }
% \author{Martin Maechler}
\seealso{
  Concrete provided copula family objects, such as
  \code{\link{copAMH}}, \code{\link{copClayton}}, etc.
}
\examples{
showClass("acopula")
str(copClayton)
}
\keyword{classes}
