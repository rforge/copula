\name{opower}
\alias{opower}
\title{Outer Power Transformation of an Archimedean Copula}
\usage{
opower(copbase, thetabase)
}
\description{ Perform the \dQuote{outer power} transformation of an
  Archimedean copula.
  %% FIXME
}
\arguments{
  \item{copbase}{a "base" copula, i.e. of class "acopula";
    must be one of the 5 five predefined families (?? FIXME ??).}
  \item{thetabase}{the (univariate) parameter 'theta' for the base
    copula.}
}
\value{a new "acopula" object; the outer power}
\author{Marius Hofert}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
% \references{
% }
\seealso{
  The class \code{\linkS4class{acopula}}; our predefined "acopula"
  family objects, \code{\link{acopula-families}}.
}
\examples{
 %% TODO
}
\keyword{manip}
%% TODO: a better \keyword{..}
