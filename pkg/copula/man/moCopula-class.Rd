\name{moCopula-class}
\docType{class}
\alias{moCopula-class}
\title{Class "moCopula" }
\description{Marshall-Olkin Copula.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("moCopula",
    ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{dimension}:}{Object of class \code{"numeric"}, the
      dimension of the copula. }
    \item{\code{exprdist}:}{Object of class \code{"expression"},
      expressions for the cdf and pdf of the copula. These expressions
      are used in function \code{pCopula()} and \code{dCopula()}.}
    \item{\code{parameters}:}{Object of class \code{"numeric"},
      parameter values. }
    \item{\code{param.names}:}{Object of class \code{"character"},
      parameter names. }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"},
      parameter lower bound. }
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"},
      parameter upper bound. }
    \item{\code{fullname}:}{Object of class \code{"character"}, family names
      of the copula (deprecated).}
  }
}
\section{Methods}{
  \describe{
    \item{dCopula}{\code{signature(copula = "moCopula")}: ... }
    \item{pCopula}{\code{signature(copula = "moCopula")}: ... }
    \item{rCopula}{\code{signature(copula = "moCopula")}: ... }
  }
}
\section{Extends}{
  Class \code{"moCopula"} extends class \code{"\linkS4class{copula}"} directly.
}
%\author{Marius Hofert}
\references{
  Nelsen, R. B. (2006), \emph{An introduction to Copulas}, Springer, New York.
}
\seealso{
  \code{\link{copula-class}}, \code{\link{moCopula-class}}.
}
\keyword{classes}
