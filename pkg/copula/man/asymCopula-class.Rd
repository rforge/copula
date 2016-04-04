\name{asymCopula-class}
\title{Class \code{"asymCopula"} and Subclasses of Asymmetric Copulas}
\docType{class}
\alias{asymCopula-class}
\alias{asymBivCopula-class}
\alias{asymExplicitCopula-class}
%
\description{
  The \emph{virtual} class \code{"asymBivCopula"} and its two (\sQuote{actual})
  subclasses \code{"asymBivCopula"} and \code{"asymExplicitCopula"}
  represent asymmetric copulas constructed from two (simple) copulas of
  the same dimension.
}
\section{Objects from the Class}{
  Objects are typically created via \code{\link{asymBivCopula}(...)}
  or \code{\link{asymExplicitCopula}(...)}.
}
\section{Slots}{
  As these classes extend \code{\linkS4class{"copula"}}, they have all
  its slots, \code{dimension}, \code{parameters}, \code{param.names},
  \code{param.lowbnd}, \code{param.upbnd}, and \code{fullname}.
  Additionally, \code{"asymCopula"} has slots
  \describe{
    \item{\code{copula1}:}{Object of class \code{"\linkS4class{copula}"}.}
    \item{\code{copula2}:}{second object of class \code{"copula"},}
  }
  and \code{"asymExplicitCopula"} has the additional slots
  \describe{
    \item{\code{exprdist}:}{an \code{\link{expression}}, ...}% FIXME
    \item{\code{derExprs1}:}{an \code{\link{expression}}, ...}% FIXME
    \item{\code{derExprs2}:}{an \code{\link{expression}}, ...}% FIXME
  }
  (whereas \code{"asymBivCopula"} has no extra slots).
}
\section{Methods}{
  Importantly, \code{\link{A}} and \code{\link{dCopula}} methods where
  devised explicitly, the latter allowing likelihood fitting (in principle).
}
%% \references{
%% }
%% \author{Ivan Kojadinovic}
\seealso{The typical \dQuote{creator} functions,
  \code{\link{asymBivCopula}()} and \code{\link{asymExplicitCopula}()}.
}
\examples{
showClass("asymCopula")# two subclasses

## all methods applicable to these subclasses:
(meths <- sapply(names(getClass("asymCopula")@subclasses),
                 function(CL) methods(class = CL),
                 simplify=FALSE))
}
\keyword{classes}
