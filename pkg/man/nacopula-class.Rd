\name{nacopula-class}
\Rdversion{1.1}
\docType{class}
\alias{outer_nacopula-class}
\alias{nacopula-class}
\alias{dim,nacopula-method}

\title{Class "nacopula" of Nested Archimedean Copulas}
\description{
  The class of Nested Archimedean Copulas, \code{"nacopula"},

  and its \emph{specialization} \code{"outer_nacopula"},
  ..........
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("nacopula", ...)},
  which is only for experts; rather they are typically constructed by

  ................
}
\section{Slots}{
  \describe{
    \item{\code{copula}:}{an object of class
      \code{"\linkS4class{acopula}"}, denoting the top-level copula of
      the nac.}
    \item{\code{comp}:}{an \code{\link{integer}} vector (possibly of
      length 0) of \bold{comp}onent indices from \code{1:d} where
      \eqn{d} is the \emph{dimension} (of the random variables), see the
      \code{dim()} method below.}
    \item{\code{childCops}:}{A \code{\link{list}} (possibly empty) of
      further nacs, i.e., objects of class \code{"nacopula"}.}
  }
}
\section{Methods}{
  \describe{
    \item{dim}{\code{signature(x = "nacopula")}: returns the dimension
      \eqn{d} of the random variables, i.e., the dimension of \eqn{u}
      when the nac describes the distribution of \eqn{U = u}.}
  }
}
% \note{ ~~further notes~~ }
\seealso{
  \code{\linkS4class{nacopula}}.
}
\examples{
showClass("nacopula")
showClass("outer_nacopula")
}
\keyword{classes}
