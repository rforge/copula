\name{psiDabsMC}
\alias{psiDabsMC}
\title{Absolute Value of Generator Derivatives via Monte Carlo}
\description{
  Computes absolute values of generator \eqn{\psi}{psi} derivatives via
  Monte Carlo simulation.
}
\usage{
psiDabsMC(t, family, theta, degree = 1, n.MC,
          method=c("log", "direct", "pois.direct", "pois"),
          log=FALSE)
}
\arguments{
  \item{t}{evaluation points}
  \item{family}{Archimedean family (name or object)}
  \item{theta}{parameter value}
  \item{degree}{order of derivative}
  \item{n.MC}{Monte Carlo sample size}
  \item{method}{different methods
    \describe{
      \item{\code{"log"}:}{proper log}
      \item{\code{"direct"}:}{direct MC formula.}
      \item{\code{"pois.direct"}:}{directly uses the Poisson density}
      \item{\code{"pois"}:}{intelligently uses the Poisson density with
	a logarithm space sum.%not exported / documented: \code{\link{lsum}()}
      }
    }
  }
  \item{log}{if TRUE the log of psiDabs is returned.}
}
\author{Marius Hofert}
\details{
... FIXME
}
\value{
... FIXME
}
\seealso{
  \code{\link{acopula-families}}.
}
\examples{
## ... FIXME
}
\keyword{distribution}
