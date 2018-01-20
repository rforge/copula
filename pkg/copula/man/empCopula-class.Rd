\name{empCopula-class}
\title{Class "empCopula" of Empirical Copulas}
\docType{class}
\alias{empCopula-class}
\description{
  Empirical Copula class.
}
\section{Objects from the Class}{
  Created by calls of the form \code{new("empCopula", ...)} or rather
  typically by \code{\link{empCopula}()} based on a matrix
  \code{X} of pseudo-observations. Implemented choices
  of smoothing are available.
}
\section{Slots}{
  \describe{
    \item{\code{X}:}{\code{\link{matrix}} of pseudo-observations based
      on which the empirical copula is constructed.}
    \item{\code{smoothing}:}{\code{\link{character}} vector determining
      the smoothing method.}
    \item{\code{offset}:}{\code{\link{numeric}} giving the shift in the
      normalizing factor for computing the empirical copula.}
    \item{\code{ties.method}:}{\code{\link{rank}()}'s ties method for
      computing the empirical copula.}
    \item{\code{dimension}:}{inherited from the
      super class \code{\linkS4class{copula}} (as selected other slots).}
  }
}
%\author{Marius Hofert}
\seealso{
  \code{\link{ellipCopula}}, \code{\link{archmCopula}},
  \code{\link{evCopula}}. And \code{\link{fhCopula}},
  \code{lowfhCopula}, \code{upfhCopula}.
}
\keyword{classes}
