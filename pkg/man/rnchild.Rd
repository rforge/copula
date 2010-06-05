\name{rnchild}
\title{Random Vector Generation For (Child) 'nacopula's}
\alias{rnchild}
\alias{rnchild-methods}
\alias{rnchild,nacopula-method}
\description{
  Methods for random vector (number) generation of nested archimedean
  copulas.
}
\usage{
rnchild(x, theta0, V0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an \R object of \code{\link{class}}
    \code{"\linkS4class{nacopula}"}.% typically produced by .....
    % once we *HAVE* a function to create these !
  }
  \item{theta0}{the parameter (vector) of the copula.}
  \item{V0}{a numeric vector of length, ...}
  \item{\dots}{further arguments, for some copula families,
    .......
    ....... %% FIXME
  }
}
\section{Methods}{
  \describe{
    \item{\code{signature(x = "nacopula")}}{ ... %% FIXME
    }
  }
}
\details{
  The generation happens recursively, descending the tree of nested
  copulas
  ... ... % FIXME
}
\value{
  a list with components
 \item{U}{a numeric matrix  .....FIXME ...}
 \item{indcol}{an index vector of column indices for \code{U}  ....FIXME ...}
%% ...
}
% \references{
% %% ~put references to the literature/web site here ~
% }
\seealso{
  \code{\link{rnacopula}},   \code{"\linkS4class{nacopula}"}. %% more
}
\examples{
## .........
}
\keyword{methods}
\keyword{distribution}

