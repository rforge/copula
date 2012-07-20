\name{setTheta}
\alias{setTheta}
\title{Specify the parameter of an Archimedean copula}
\description{
  Set or change the parameter \eqn{\theta}{theta} (\code{theta}) of
  an Archimedean copula \code{x} of class \code{"\linkS4class{acopula}"}.  This
  is used for constructing copula models with specified parameter, as, for
  example, in \code{\link{onacopula}()}.
}
\usage{
setTheta(x, value, na.ok = TRUE, noCheck = FALSE)
}
\arguments{
  \item{x}{acopula}
  \item{value}{parameter value, \code{\link{numeric}} or
    \code{\link{NA}} (when \code{na.ok} is true.)}
  \item{na.ok}{logical indicating if \code{\link{NA}} values are ok for
    \code{theta}.}
  \item{noCheck}{logical indicating if parameter contraint checks should
    be skipped.}
}
\value{
  an \code{"\linkS4class{acopula}"} as \code{x}, just with the slot
  \code{theta} set to \code{value}.
}
\author{Martin Maechler}
\examples{
myC <- setTheta(copClayton, 0.5)
myC

# negative theta is ok for dim = 2 :
myF <- setTheta(copFrank, -2.5, noCheck=TRUE)
myF@tau(myF@theta) # -0.262
}
\keyword{manip}
