\name{ACopula-class}
\Rdversion{1.1}
\docType{class}
\alias{ACopula-class}

\title{Class "ACopula" of Archimedean Copula Families}
\description{
  This class \code{"ACopula"} Archimedean Copula Families
  is mainly used for providing objects of known A.copula families
  with all related functions.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ACopula", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{a string (class \code{"character"}) describing
      the copula family.}
    \item{\code{psi}, \code{psiInv}:}{\eqn{psi}, the generator
      \code{\link{function}} and its inverse.}
    ## to do
    \item{\code{paraConstr}:}{Object of class \code{"function"} ~~ }
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
No methods defined with class "ACopula" in the signature.
}
% \references{ ~put references to the literature/web site here ~ }
% \author{Martin Maechler}
\seealso{
  Concrete provided copula family objects, such as
  \code{\link{copAMH}}, \code{\link{copClayton}}, etc.
}
\examples{
showClass("ACopula")
}
\keyword{classes}
