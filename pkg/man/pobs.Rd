\name{pobs}
\alias{pobs}
\title{Pseudo-observations}
\description{
  Computes the pseudo-observations for the given data matrix
}
\usage{
pobs(x, na.last = "keep",
     ties.method = c("average", "first", "random", "max", "min"))
}
\arguments{
  \item{x}{numeric matrix of random variates to be converted to
    pseudo-observations.}
  \item{na.last, ties.method}{are passed to \code{\link{rank}()}, see
    there.}
}
\value{matrix of the same dimensions as \code{x} with the
  pseudo-observations.
}
\author{Marius Hofert}
\examples{
## Simple definition of the function:
pobs

## TODO
}


