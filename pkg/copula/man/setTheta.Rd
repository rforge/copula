\name{setTheta}
\alias{setTheta}
\alias{setTheta,acopula,ANY-method}
\alias{setTheta,outer_nacopula,numeric-method}
\title{Specify the Parameter of a (Nested) Archimedean Copula}
\description{
%% FIXME:
  Set or change the parameter \eqn{\theta}{theta} (\code{theta}) of
  a (nested) Archimedean copula \code{x} of class
  \code{"\linkS4class{acopula}"} or
  \code{"\linkS4class{outer_nacopula}"}.  This is used for constructing
  copula models with specified parameter, as, for example, in
  \code{\link{onacopula}()}.
}
\usage{
setTheta(x, value, \dots)
\S4method{setTheta}{acopula,ANY}(x, value, na.ok = TRUE, noCheck = FALSE)
\S4method{setTheta}{outer_nacopula,numeric}(x, value, na.ok = TRUE, noCheck = FALSE)
}
\arguments{
  \item{x}{acopula}
  \item{value}{parameter value, \code{\link{numeric}} or
    \code{\link{NA}} (when \code{na.ok} is true.)}
  \item{\dots}{further arguments for methods.}
  \item{na.ok}{logical indicating if \code{\link{NA}} values are ok for
    \code{theta}.}
  \item{noCheck}{logical indicating if parameter contraint checks should
    be skipped.}
}
\value{
  an \R object of the same class as \code{x}, with the main parameter
  (often called \code{theta}) set to \code{value}.
%% FIXME: be more explicit about outer_nacopula  case
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
