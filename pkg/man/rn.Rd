\name{rn}
\alias{rn}
\alias{rn-methods}
\alias{rn,outer_nACopula-method}
\title{Random Vector Generation For 'outer_nACopula's}
\description{
  Methods for random vector (number) generation of nested archimedean
  copulas.
}
\usage{
rn(x, n, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an \R object of \code{\link{class}}
    \code{"\linkS4class{outer_nACopula}"}.% typically produced by .....
    % once we *HAVE* a function to create these !
  }
  \item{n}{integer specifying the sample size, i.e., the number of
    copula-distributed random vectors \eqn{\mathbf{U}_i}, to be generated.
  }
  \item{\dots}{further arguments, for some copula families,
    .......
    ....... %% FIXME
  }
}
\details{
  The generation happens by calling \code{\link{rnchild}()} on
  each child copula (which itself recursively descends the tree of nested
  copulas)
  ... ... % FIXME
}
\value{
  a numeric matrix, say \code{U}, ... .. % FIXME
}
\section{Methods}{
  \describe{
    \item{\code{signature(x = "outer_nACopula")}}{
    }
  }
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
% \references{
% %% ~put references to the literature/web site here ~
% }
\seealso{
  \code{\link{rnchild}}; classes \code{"outer_nACopula"} and
  \code{"\linkS4class{nACopula}"}.
}
\examples{
}
\keyword{methods}
\keyword{distribution}

