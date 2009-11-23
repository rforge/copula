\name{nACopula-class}
\Rdversion{1.1}
\docType{class}
\alias{final_nACopula-class}
\alias{nACopula-class}
\alias{dim,nACopula-method}

\title{Class "nACopula" of Nested Archimedean Copulas}
\description{
  The class of Nested Archimedean Copulas, \code{"nACopula"},

  and its \emph{specialization} \code{"final_nACopula"},
  ..........
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("nACopula", ...)},
  which is only for experts; rather they are typically constructed by

  ................
}
\section{Slots}{
  \describe{
    \item{\code{copula}:}{an object of class
      \code{"\linkS4class{ACopula}"}, denoting the top-level copula of
      the NAC.}
    \item{\code{comp}:}{an \code{\link{integer}} vector (possibly of
      length 0) of \bold{comp}onent indices from \code{1:d} where
      \eqn{d} is the \emph{dimension} (of the random variables), see the
      \code{dim()} method below.}
    \item{\code{childCops}:}{A \code{\link{list}} (possibly empty) of
      further NACs, i.e., objects of class \code{"nACopula"}.}
  }
}
\section{Methods}{
  \describe{
    \item{dim}{\code{signature(x = "nACopula")}: returns the dimension
      \eqn{d} of the random variables, i.e., the dimension of \eqn{u}
      when the NAC describes the distribution of \eqn{U = u}.}
  }
}
% \note{ ~~further notes~~ }
\seealso{
  \code{\linkS4class{ACopula}}.
}
\examples{
showClass("nACopula")
showClass("final_nACopula")
}
\keyword{classes}
