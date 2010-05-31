\name{prob}
\alias{prob}
\alias{prob-methods}
\alias{prob,outer_nacopula-method}
\title{Probability of "nacopula" (in a Hypercube)}
\description{
  Compute probabilities \eqn{U \in [l,u]}{U in [l,u]}, where U is
  (\eqn{d}-dimensionally) distributed according to the nested
  archimedean copula \code{x}, and \eqn{l} and \eqn{u} are
  \eqn{d-}dimensional vectors of \bold{l}ower and \bold{u}pper interval
  boundaries, together defining a \eqn{d}-dimensional hypercube.
}
\usage{
prob(x, l, u)
}
\arguments{
  \item{x}{a (nested) Archimedean copula object, i.e., from class
    \code{\linkS4class{nacopula}}, typically the result of \code{\link{onacopula}(..)}.}
  \item{l, u}{numeric \eqn{d}-dimensional interval boundaries,
    \eqn{0 \le l_i \le u_i \le 1}{0 <= l_i <= u_i <= 1}, for \eqn{i = 1,\dots,d}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\section{Methods}{
  \describe{
    \item{\code{signature(x = "outer_nacopula")}}{
    }
  }
}
\value{a number (in \eqn{[0,1]}), the probability
  \eqn{Pr[l_i \le U_i \le u_i]}{Pr[l[i] <= U[i] <= u[i]]}.
}
% \references{
% %% ~put references to the literature/web site here ~
% }
% \author{Marius Hofert}
\seealso{
  \code{\link{onacopula}()} for construction;
  the class, \code{\linkS4class{nacopula}}, and e.g., evaluation of
  copulas, \code{\link{pnacopula}(.)}.
}
\examples{
C3d.AMH <- onacopula("AMH", C(0.7135, 1, C(0.943, 2:3)))
prob(C3d.AMH, l = c(.1, .2, .3),
              u = c(.3, .3, .5))
##-> 0.005662
}
\keyword{ ~kwd1 }
\keyword{methods}

