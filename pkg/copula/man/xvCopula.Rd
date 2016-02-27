\name{xvCopula}
\title{Model (copula) selection based on \code{k}-fold cross-validation}
\alias{xvCopula}
\description{
  Computes the leave-one-out cross-validation criterion (or a
  \code{k}-fold version of it) for the hypothesized parametric copula
  family based on maximum pseudo-likelihood estimation. The
  leave-one-out criterion is given in equation (42) of
  Gr\enc{ø}{o}nneberg and Hjort (2014). Its \code{k}-fold version
  results from considering \code{k} data blocks instead of \code{n},
  where \code{n}
  is the sample size. When \code{n} is very large, the latter can be
  significantly faster if \code{k} is "small" compared to \code{n}.
}
\usage{
xvCopula(copula, x, k = NULL, verbose = TRUE, \dots)
}
\arguments{
  \item{copula}{object of class \code{"\linkS4class{copula}"} representing the
    hypothesized copula family.}
  \item{x}{a data matrix that will be transformed to pseudo-observations.}
  \item{k}{the number of data blocks; if \code{k = NULL}, \code{nrow(x)}
    blocks are considered (which corresponds to leave-one-out
    cross-validation).}
  \item{\dots}{additional arguments passed to \code{fitCopula()}.}
  }
%\details{}
\value{
  A real number equal to the cross-validation criterion
  multiplied by the sample size.
}
\references{
  Gr\enc{ø}{o}nneberg, S., and Hjort, N.L. (2014), The copula information
  criteria, \emph{Scandinavian Journal of Statistics}
  \bold{41}, 436--459.
}
%\note{}
\seealso{
  \code{\link{fitCopula}()} for the underlying estimation procedure and
  \code{\link{gofCopula}()} for goodness-of-fit tests.
}
\examples{
\donttest{
## A two-dimensional data example ----------------------------------
x <- rCopula(200, claytonCopula(3))

## Model (copula) selection
xvCopula(gumbelCopula(), x)
xvCopula(frankCopula(), x)
xvCopula(joeCopula(), x)
xvCopula(claytonCopula(), x)
xvCopula(normalCopula(), x)
xvCopula(tCopula(), x)
xvCopula(plackettCopula(), x)


## The same with k-fold cross-validation
k <- 50
xvCopula(gumbelCopula(), x, k=k)
xvCopula(frankCopula(), x, k=k)
xvCopula(joeCopula(), x, k=k)
xvCopula(claytonCopula(), x, k=k)
xvCopula(normalCopula(), x, k=k)
xvCopula(tCopula(), x, k=k)
xvCopula(plackettCopula(), x, k=k)


}% dont..
}
\keyword{models}
\keyword{multivariate}
