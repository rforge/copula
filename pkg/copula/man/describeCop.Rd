\name{describeCop}
\title{Copula Short Description as String}
\alias{describeCop}
\alias{describeCop-methods}
\alias{describeCop,xcopula,ANY-method}
\alias{describeCop,copula,missing-method}
\alias{describeCop,archmCopula,character-method}
\alias{describeCop,copula,character-method}
\alias{describeCop,ellipCopula,character-method}
\alias{describeCop,fgmCopula,character-method}
\docType{methods}
\description{
  Describe a \code{\linkS4class{copula}} object, i.e., its basic
  properties as a string.  This is a utility used when
  \code{\link{print}()}ing or \code{\link{plot}()}ting copulas, e.g.,
  after a fitting.
}
\usage{
describeCop(x, kind = c("short", "very short", "long"))
}
\section{Methods}{
  \describe{
    \item{\code{signature(x = "archmCopula", kind = "ANY")}}{ .. }
    \item{\code{signature(x = "copula", kind = "character")}}{ .. }
    \item{\code{signature(x = "copula", kind = "missing")}}{ .. }
    \item{\code{signature(x = "ellipCopula", kind = "character")}}{ .. }
    \item{\code{signature(x = "fgmCopula", kind = "ANY")}}{ .. }
    \item{\code{signature(x = "xcopula", kind = "ANY")}}{ .. }
  }
}
\arguments{
  \item{x}{a \code{\linkS4class{copula}} object, or a generalization
    such as  \code{\linkS4class{parCopula}}.}
  \item{kind}{a \code{\link{character}} string specifying the size (or
    \dQuote{complexity} of the copula description desired.}
}
%% \details{
%% }
\value{
  a \code{\link{character}} string.
}
\seealso{
  Copula class definition \code{\linkS4class{copula}};
}
\examples{
## FIXME
}
\keyword{utilities}
\keyword{methods}
