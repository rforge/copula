\name{interval-class}
\Rdversion{1.1}
\docType{class}
\alias{interval-class}
\alias{maybeInterval-class}
%
\alias{format,interval-method}
\alias{show,interval-method}
\alias{Summary,interval-method}

\title{Class "interval" of simple Intervals}
\description{
  The S4 \code{\link{class}}) \code{"interval"} is a simple class for numeric
  interval.

  \code{"maybeInterval"} is a simple class union (see
  \code{\link{setClassUnion}}) of \code{"interval"} and \code{"NULL"}.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("interval", ...)},
  but typically rather via \code{\link{interval}()}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{numeric vector of length two, specifying the
      interval ranges.}
    \item{\code{open}:}{\code{\link{logical}} vector of length two,
      specifying if the interval is open or closed on the left and
      right, respectively.}
  }
}
\section{Extends}{
  Class \code{"interval"} extends \code{"\linkS4class{numeric}"}, from
  data part, and \code{"maybeInterval"}, directly.
}
\section{Methods}{
  \describe{
    \item{format}{\code{signature(x = "interval")}: ... }
    \item{show}{\code{signature(object = "interval")}: ... }
    \item{Summary}{\code{signature(x = "interval")}: Group methods,
      notably \code{\link{range}()}, \code{\link{min}()}, etc.}
  }
}
\author{Martin Maechler}
\note{
  There are more sophisticated interval classes, functions and methods,
  notably in package \pkg{Intervals}.  We only use this as a simple
  interface in order to specify our copula functions consistently.
}
\seealso{
  \code{\link{interval}} constructs "interval" objects conveniently.
}
\keyword{classes}
