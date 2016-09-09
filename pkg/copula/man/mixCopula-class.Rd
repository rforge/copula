\name{mixCopula-class}
\title{Class \code{"mixCopula"} of Copula Mixtures}
\alias{mixCopula-class}
\alias{dim,mixCopula-method}
\alias{lambda,mixCopula-method}
\alias{rho,mixCopula-method}
\Rdversion{1.1}
\docType{class}
\description{
  The class \code{"mixCopula"} is the class of all finite mixtures of
  copulas.

  These are given by (a list of) \eqn{m} \dQuote{arbitrary} copulas, and
  their respective \eqn{m} non-negative weights.
}
\section{Objects from the Class}{
  Objects are typically created by \code{\link{mixCopula}()}.
}
\section{Slots}{
  \describe{
    \item{\code{w}:}{Object of class \code{"mixWeights"} ~~ }
    \item{\code{cops}:}{Object of class \code{"parClist"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{parCopula}"}, directly.
Class \code{"\linkS4class{Copula}"}, by class "parCopula", distance 2.
}
\section{Methods}{
  \describe{
    \item{dim}{\code{signature(x = "mixCopula")}: dimension of copula.}
    \item{rho}{\code{signature(x = "mixCopula")}: Spearman's rho of
      copula \code{x}.}
    \item{lambda}{\code{signature(x = "mixCopula")}: lower and upper
      tail dependecies \code{\link{lambda}}, \eqn{(\lambda[L],\lambda[U])},
      of the mixture copula.}
  }
}
\seealso{
  \code{\link{mixCopula}} for creation and examples.
}
\examples{
showClass("mixCopula")
}
\keyword{classes}
