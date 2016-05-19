\name{fixParam}
\title{Fix a Subset of a Copula Parameter Vector}
\alias{fixParam}
\alias{fixedParam<-}
\description{
  It is sometimes useful to keep fixed some components of a copula
  parameter vector whereas the others are \dQuote{free} and will be
  estimated, e.g., by \code{\link{fitCopula}}.
}
\usage{
fixParam(param, fixed = TRUE)
fixedParam(copula) <- value
}
\arguments{
  \item{param}{numeric parameter vector}
  \item{fixed}{logical vector of the same length as param: TRUE = fixed}
  \item{copula}{a copula object}
  \item{value}{a logical vector to modify the fixedness of the parameters}
}
\value{
  A numerical vector with attribute \code{"fixed"} (a
  \code{\link{logical}}, either \code{TRUE} or vector of the same length
  as \code{param}) to indicate which components of \code{param} are to
  be held fixed or not.
}
%\author{Jun Yan}
%% \details{
%% }
%% \references{
%% }
\seealso{
  \code{\link{fitCopula}} for fitting;
  %% FIXME  t-copulas (\code{\link{tCopula}}) ...
}
\examples{
normalCopula(dim = 3, fixParam(c(.6,.3,.2), c(TRUE, FALSE,FALSE)),
             dispstr = "un")
%% FIXME: print() / show() should indicated "fixed"
}
\keyword{manip}
\keyword{utilities}
