\name{Cn}
\alias{Cn}
\title{The Empirical Copula}
\description{
  Given a random sample from a distribution with continuous margins and
  copula C, the empirical copula is a natural nonparametric estimator of
  C. The function \code{Cn} computes the empirical copula.
}
\usage{
Cn(u, U, do.pobs=TRUE, offset=0, method=c("C", "R"))
}
\arguments{
  \item{u}{a \code{\link{matrix}} with elements in \eqn{[0,1]} whose
    lines contain the evaluation points of the empirical copula.}
  \item{U}{a \code{\link{matrix}} with elements in \eqn{[0,1]} and
    with the same number of colums as \code{u}. The lines of \code{U}
    contain (pseudo-)data based on which the empirical copula is built.
    A multivariate random sample can be transformed to an appropriate
    \code{U} via \code{\link{pobs}()}. This can also be achieved by
    \code{do.pobs=TRUE}.}
  \item{do.pobs}{\code{\link{logical}} indicating whether
    \code{\link{pobs}()} is applied to \code{U} before the empirical
    copula is computed (defaults to \code{TRUE}).}
  \item{offset}{scaling factor (defaults to 0) of the form \code{sum(..)/(n+offset)}
    where \code{n} denotes the number of rows of \code{U}.}
  \item{method}{method string, indicating which method is applied to
    compute the empirical copula.}
}
\value{
  Returns a vector of length equal to the number of rows of \code{u}
  containing the values of the empirical copula of \code{U} at \code{u}.
}
\details{
  There are several asymptotically equivalent definitions of the
  empirical copula. Here, the empirical copula is simply defined as
  the empirical distribution function computed from the
  pseudo-observations, that is,
  \deqn{C(\bm{u})=\frac{1}{n+1}\sum_{k=1}^n\mathbf{1}_{\{\hat{\bm{U}}_i\le\bm{u}\}}}{1/(n+1)sum(I_{U_i<=u},
  i=1,..,n)},
  where \eqn{\hat{\bm{U}}_i}{U_i}, \eqn{i\in\{1,\dots,n\}}{i=1,..,n}, denote the
  pseudo-observations (rows in \code{U}) and \eqn{n} the sample size
  (number of rows of \code{U}).
}
\author{Ivan Kojadinovic, Marius Hofert}
\references{
  R\enc{ü}{u}schendorf, L. (1976).  Asymptotic distributions of
  multivariate rank order statistics, \emph{Annals of
    Statistics}. \bold{4}, 912--923.

  Deheuvels, P. (1979).  La fonction de dépendance empirique et ses
  propriétés: un test non paramétrique d'indépendance, \emph{Acad. Roy.
    Belg. Bull. Cl. Sci.}, 5th Ser. \bold{65}, 274--292.

  Deheuvels, P. (1981).  A non parametric test for independence,
  \emph{Publ. Inst. Statist. Univ. Paris}. \bold{26}, 29--50.
}
\seealso{\code{\link{pobs}()} for computing pseudo-observations,
  \code{\link{pCopula}()} for evaluating a copula.}
\examples{
n <- 100
d <- 3
cop <- onacopulaL("Gumbel", list(theta=2, 1:d))
set.seed(1)
U <- rCopula(n, cop)

## random points were to evaluate the empirical copula
u <- matrix(runif(n*d), n, d)
ec <- Cn(u, U=U)

## compare with true distribution function
mean(abs(pCopula(u, copula=cop)-ec)) # increase n to decrease this error

## compare the empirical copula and the true copula
## on the diagonal of the unit square
Cn. <- function(x) Cn(do.call(cbind, rep(list(x), d)), U=U)
curve(Cn., 0, 1, main="Diagonal of a Gumbel copula",
      xlab="u", ylab=expression(italic(C)[n](italic(u),..,italic(u))))
pC <- function(x) pCopula(do.call(cbind, rep(list(x), d)), cop)
curve(pC, lty=2, add=TRUE)
legend("topleft", lty=1:2, bty="n",inset=0.02,
       legend=c(expression(italic(C)[n]), expression(italic(C))))

## check the empirical copula with its Kendall distribution function
plot( pK(Cn(U, U=U), cop=cop@copula, d=d) ) # must be uniform
}
\keyword{empirical copula}
