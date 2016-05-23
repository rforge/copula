\name{khoudrajiCopula-class}
\title{Class \code{"khoudrajiCopula"} and its Subclasses}
\docType{class}
\alias{khoudrajiCopula-class}
\alias{khoudrajiBivCopula-class}
\alias{khoudrajiExplicitCopula-class}

\description{
  The \emph{virtual} class \code{"khoudrajiCopula"} and its two (\sQuote{actual})
  subclasses \code{"khoudrajiBivCopula"} and \code{"khoudrajiExplicitCopula"}
  represent copulas constructed using Khoudraji's device from two copulas of
  the same dimension.
}
\section{Objects from the Class}{
  Objects are typically created via \code{\link{khoudrajiCopula}(...)}.
}
\section{Slots}{
  As these classes extend \code{"\linkS4class{copula}"}, they have all
  its slots, \code{dimension}, \code{parameters}, \code{param.names},
  \code{param.lowbnd}, \code{param.upbnd}, and \code{fullname}.
  The classes \code{"khoudrajiCopula"} and \code{"khoudrajiBivCopula"}
  have the extra slots
  \describe{
    \item{\code{copula1}:}{Object of class \code{"\linkS4class{copula}"}.}
    \item{\code{copula2}:}{second object of class \code{"copula"},}
  }
  while the class \code{"khoudrajiExplicitCopula"} has the additional slots
  \describe{
    \item{\code{exprdist}:}{an \code{\link{expression}}, ...}% FIXME
    \item{\code{derExprs1}:}{an \code{\link{expression}} of length \eqn{d}, ...}% FIXME
    \item{\code{derExprs2}:}{an \code{\link{expression}} of length \eqn{d}, ...}% FIXME
  }
}
\section{Methods}{
  TODO:%Importantly, \code{\link{A}} and \code{\link{dCopula}} methods where
  %devised explicitly, the latter allowing likelihood fitting (in principle).
}
%% \references{
%% }
%% \author{Ivan Kojadinovic}
\seealso{ \code{\link{khoudrajiCopula}()}}
\examples{
showClass("khoudrajiCopula")# two subclasses

## all methods applicable to these subclasses:
(meths <- sapply(names(getClass("khoudrajiCopula")@subclasses),
                 function(CL) methods(class = CL),
                 simplify=FALSE))
}
\keyword{classes}
