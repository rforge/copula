\name{dnacopula}
\alias{dnacopula}
\title{Copula Density Evaluation}
\description{
  Returns the copula density at u
}
\usage{
dnacopula(x, u, MC, log=FALSE)
}
\arguments{
  \item{x}{an object of class \code{"\linkS4class{outer_nacopula}"}.}
  \item{u}{argument of the copula x}
  \item{MC}{(optional) integer or \code{NULL}.  If specified (not NULL)
    Monte Carlo is used for evaluation of the density with sample size
    \eqn{=} \code{MC}.}
  \item{log}{if TRUE the log-density is evaluated}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  numeric ...
}
% \references{
% %% ~put references to the literature/web site here ~
% }
\seealso{
 ...
}
\examples{
}
\keyword{distribution}

