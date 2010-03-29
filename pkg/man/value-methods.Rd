\name{value-methods}
\docType{methods}
\alias{value}
\alias{value-methods}
\alias{value,nACopula-method}
\title{Nested A.Copula Evaluation -- value() Generic Function and Methods}
\description{
  \code{value(<NAC>, u)} evaluates a nested Archimedean copula (an
  object of class \code{\linkS4class{nACopula}}) at \code{u}.
}
\usage{
value(x, u)
}
\arguments{
  \item{x}{an \R object \dQuote{at which evaluation should happen};
    currently must be of class
    \code{\linkS4class{nACopula}}, typically of (subclass)
    \code{\linkS4class{outer_nACopula}}.}
  \item{u}{numeric vector length \eqn{d}.}
}
\details{
.... %%
}
\value{
  a number (for now; if parallelized, allowing \code{u} to be a
  \eqn{n\times d}{n * d} matrix, the return value will be a numeric
  (vector) of length \eqn{n}).
}
\references{
....%% ~put references to the literature/web site here ~
}
\section{Methods}{
  \describe{
    \item{\code{signature(x = "nACopula")}}{ ... (see above) ... }
  }
}
\examples{
%% TODO
}
\keyword{methods}

