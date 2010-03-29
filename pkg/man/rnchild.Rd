\name{rnchild}
\title{Random Vector Generation For (Child) 'nACopula's}
\alias{rnchild}
\alias{rnchild-methods}
\alias{rnchild,nACopula-method}
\description{
  Methods for random vector (number) generation of nested archimedean
  copulas.
}
\usage{
rnchild(x, n, psi0Inv, theta0, V0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an \R object of \code{\link{class}}
    \code{"\linkS4class{nACopula}"}.% typically produced by .....
    % once we *HAVE* a function to create these !
  }
  \item{n}{integer specifying the sample size, i.e., the number of
    copula-distributed random vectors \eqn{\mathbf{U}_i}, to be generated.
  }
  \item{psi0Inv}{
    the function \eqn{\psi_0^{-1}}{psi0^(-1)} for ... ... %FIXME
  }
  \item{theta0}{the parameter (vector) of the copula.}
  \item{V0}{a numeric vector of length \code{n}, ...}
  \item{\dots}{further arguments, for some copula families,
    .......
    ....... %% FIXME
  }
}
\section{Methods}{
  \describe{
    \item{\code{signature(x = "nACopula")}}{ ... %% FIXME
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
  \code{\link{rn}},   \code{"\linkS4class{nACopula}"}. %% more
}
\examples{
## .........
}
\keyword{methods}
\keyword{distribution}

