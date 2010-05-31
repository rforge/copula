\name{pnacopula-methods}
\docType{methods}
\alias{pnacopula}
\alias{pnacopula-methods}
\alias{pnacopula,nacopula-method}
\title{Nested A.Copula Evaluation -- pnacopula() Generic Function and Methods}
\description{
  \code{pnacopula(<nac>, u)} evaluates a nested Archimedean copula (an
  object of class \code{\linkS4class{nacopula}}) at \code{u}.
}
\usage{
pnacopula(x, u)
}
\arguments{
  \item{x}{an \R object \dQuote{at which evaluation should happen};
    currently must be of class
    \code{\linkS4class{nacopula}}, typically of (subclass)
    \code{\linkS4class{outer_nacopula}}.}
  \item{u}{numeric vector length \eqn{d}.}
}
\details{
.... %% FIXME
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
    \item{\code{signature(x = "nacopula")}}{ ... (see above) ... }
  }
}
\examples{
%% TODO
}
\keyword{methods}

