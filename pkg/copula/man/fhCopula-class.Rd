\name{fhCopula-class}
\title{Class "fhCopula" of Fr\enc{é}{e}chet-Hoeffding Bound Copulas}
\docType{class}
\alias{fhCopula-class}
\alias{lowfhCopula-class}
\alias{upfhCopula-class}
\description{
  Fr\enc{é}{e}chet-Hoeffding bound copula class.
}
\section{Objects from the Class}{
  Created by calls of the form \code{new("fhCopula", ...)} or rather
  typically by \code{\link{fhCopula}()}. Implemented
  ``families'' are the lower and upper Fr\enc{é}{e}chet-Hoeffding bound
  copulas (\code{lowfhCopula()}, \code{upfhCopula()}).
}
\section{Slots}{
  \describe{
    \item{\code{dimension}:}{inherited from the
      super class \code{\linkS4class{copula}} (as other slots).}
    \item{\code{exprdist}:}{Object of class \code{"expression"}:
      expressions of the cdf (nonexisting pdf is \code{expression()})
      of the copula.}
  }
}
%\author{Marius Hofert}
\seealso{
  \code{\link{ellipCopula}}, \code{\link{archmCopula}},
  \code{\link{evCopula}}. And \code{\link{fhCopula}},
  \code{lowfhCopula}, \code{upfhCopula}.
}
\keyword{classes}
