\name{archmCopula-class}
\docType{class}
\alias{archmCopula-class}
\alias{claytonCopula-class}
\alias{frankCopula-class}
\alias{gumbelCopula-class}
\alias{amhCopula-class}
\alias{joeCopula-class}
\title{Class "archmCopula" }
\description{Archimedean copula class.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("archmCopula",
    ...)} or by function \code{\link{archmCopula}}.  Implemented
  families are Clayton, Gumbel, Frank, Joe, and Ali-Mikhail-Haq.
}
\section{Slots}{
  \describe{
    \item{\code{exprdist}:}{Object of class \code{"expression"}:
      expressions of the cdf and pdf of the copula. These expressions
      are used in function \code{\link{pcopula}} and \code{\link{dcopula}}. }
    \item{\code{dimension, parameters}, etc:}{all inherited from the
      super class \code{\linkS4class{copula}}.}
  }
}
\section{Methods}{
  \describe{
    \item{dcopula}{\code{signature(copula = "claytonCopula")}: ... }
    \item{pcopula}{\code{signature(copula = "claytonCopula")}: ... }
    \item{rcopula}{\code{signature(copula = "claytonCopula")}: ... }
    \item{dcopula}{\code{signature(copula = "frankCopula")}: ... }
    \item{pcopula}{\code{signature(copula = "frankCopula")}: ... }
    \item{rcopula}{\code{signature(copula = "frankCopula")}: ... }
    \item{dcopula}{\code{signature(copula = "gumbelCopula")}: ... }
    \item{pcopula}{\code{signature(copula = "gumbelCopula")}: ... }
    \item{rcopula}{\code{signature(copula = "gumbelCopula")}: ... }
    \item{dcopula}{\code{signature(copula = "amhCopula")}: ... }
    \item{pcopula}{\code{signature(copula = "amhCopula")}: ... }
    \item{rcopula}{\code{signature(copula = "amhCopula")}: ... }
    \item{dcopula}{\code{signature(copula = "joeCopula")}: ... }
    \item{pcopula}{\code{signature(copula = "joeCopula")}: ... }
    \item{rcopula}{\code{signature(copula = "joeCopula")}: ... }
  }
}
\section{Extends}{
  Class \code{"archmCopula"} extends class \code{"\linkS4class{copula}"}
  directly.  Class \code{"claytonCopula"}, \code{"frankCopula"},
  \code{"gumbelCopula"}, \code{"amhCopula"} and \code{"joeCopula"}
  extends class \code{"archmCopula"} directly.
}
\note{
  \code{"gumbelCopula"} is also of class \code{"\linkS4class{evCopula}"}.
}
\seealso{
  \code{\link{archmCopula}}, for constructing such copula objects;
  \code{\link{copula-class}}.
}
\keyword{classes}
