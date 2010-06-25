\name{getAcop}
\alias{getAcop}
\title{Get "acopula" Family Object by Name}
\description{
  Get one of our "acopula" family objects (see \link{acopula-families} by name.
}
\usage{
getAcop(family, check=TRUE)
}
\arguments{
  \item{family}{character string: short or longer form of copula family
    name.}
  \item{check}{logical indicating if the class of the return value
    should be checked to be \code{"\linkS4class{acopula}"}.}
}
\usage{
getAcop(family, check = TRUE)
}
\value{one of our \code{"\linkS4class{acopula}"} objects}
\author{Martin Maechler}
\seealso{
  Our predefined \code{\link{acopula-families}};
  the class definition \code{"\linkS4class{acopula}"}.
}
\examples{
getAcop("Gumbel")

stopifnot(identical(getAcop("J"),
                    getAcop("Joe")))
}
\keyword{manip}
