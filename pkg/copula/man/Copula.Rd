\name{Copula}
\alias{Copula}
\alias{dcopula}
\alias{pcopula}
\alias{rcopula}
% dcopula
\alias{dcopula,claytonCopula-method}
\alias{dcopula,frankCopula-method}
\alias{dcopula,gumbelCopula-method}
\alias{dcopula,normalCopula-method}
\alias{dcopula,tCopula-method}
\alias{dcopula,galambosCopula-method}
\alias{dcopula,amhCopula-method}
\alias{dcopula,joeCopula-method}
\alias{dcopula,huslerReissCopula-method}
\alias{dcopula,tawnCopula-method}
\alias{dcopula,tevCopula-method}
\alias{dcopula,plackettCopula-method}
\alias{dcopula,fgmCopula-method}
\alias{dcopula,indepCopula-method}
\alias{dcopula,asymCopula-method}
\alias{dcopula,asymExplicitCopula-method}
\alias{dcopula,nacopula-method}
% pcopula
\alias{pcopula,claytonCopula-method}
\alias{pcopula,frankCopula-method}
\alias{pcopula,gumbelCopula-method}
\alias{pcopula,normalCopula-method}
\alias{pcopula,tCopula-method}
\alias{pcopula,galambosCopula-method}
\alias{pcopula,amhCopula-method}
\alias{pcopula,joeCopula-method}
\alias{pcopula,huslerReissCopula-method}
\alias{pcopula,tawnCopula-method}
\alias{pcopula,tevCopula-method}
\alias{pcopula,plackettCopula-method}
\alias{pcopula,fgmCopula-method}
\alias{pcopula,indepCopula-method}
\alias{pcopula,asymCopula-method}
\alias{pcopula,asymExplicitCopula-method}
\alias{pcopula,nacopula-method}
% rcopula
\alias{rcopula,claytonCopula-method}
\alias{rcopula,frankCopula-method}
\alias{rcopula,gumbelCopula-method}
\alias{rcopula,normalCopula-method}
\alias{rcopula,tCopula-method}
\alias{rcopula,galambosCopula-method}
\alias{rcopula,amhCopula-method}
\alias{rcopula,joeCopula-method}
\alias{rcopula,huslerReissCopula-method}
\alias{rcopula,plackettCopula-method}
\alias{rcopula,fgmCopula-method}
\alias{rcopula,evCopula-method}
\alias{rcopula,indepCopula-method}
\alias{rcopula,asymCopula-method}
\alias{rcopula,asymExplicitCopula-method}
\alias{rcopula,nacopula-method}
% main
\title{Copula functions}
\description{
  Density, distribution function, and random generation for a
  \code{copula} object.
}
\usage{
dcopula(copula, u, log=FALSE, \dots)
pcopula(copula, u, \dots)
rcopula(copula, n, \dots)
}
\arguments{
  \item{copula}{an \R object of class \code{"\linkS4class{Copula}"},
    (i.e., \code{"\linkS4class{copula}"} or  \code{"\linkS4class{nacopula}"}).}
  \item{u}{a vector of the copula dimension or a matrix with number of
    columns being equal to the copula dimension, giving the coordinates of the
    points where the density or distribution function needs to be evaluated.
  }
  \item{log}{logical indicating if the \eqn{\log(f(\cdot))}{log(f(.))}
    should be returned instead of \eqn{f(\cdot)}{f(.)}.}
  \item{n}{(for \code{rcopula()}:) number of observations to be generated.}
  \item{\dots}{further optional arguments for some methods, e.g.,
    \code{method}.}
}
\details{
  The density (\code{dcopula}) and distribution function
  (\code{pcopula}) methods for Archimedean copulas now use the
  corresponding function slots of the Archimedean copula objects, such
  as \code{\link{copClayton}}, \code{\link{copGumbel}}, etc.

  The distribution function of a \code{t} copula uses \code{pmvt} from
  package \pkg{mvtnorm}; similarly, the density (\code{dcopula}) calls
  \code{\link[mvtnorm]{dmvt}} from \pkg{mvtnorm}.
  The \code{normalCopula} methods use \code{dmvnorm} and \code{pmvnorm}
  from the same package.

  The random number generator for an Archimedean copula uses the
  conditional approach for the bivariate case and the Marshal-Olkin (1988)
  approach for dimension greater than 2.
}
\value{
  \code{"dcopula"} gives the density, \code{"pcopula"} gives the
  distribution function, and \code{"rcopula"} generates random variates.
}
\references{
  Frees, E. W.  and Valdez, E. A. (1998). Understanding relationships using
  copulas.
  \emph{North American Actuarial Journal} \bold{2}, 1--25.

  Genest, C. and Favre, A.-C. (2007).
  Everything you always wanted to know about copula modeling but were
  afraid to ask.
  \emph{Journal of Hydrologic Engineering} \bold{12}, 347--368.

  Joe, H. (1997). \emph{Multivariate Models and Dependence Concepts}.
  Chapman and Hall, London.

  Marshall, A. W. and Olkin, I. (1988)
  Families of multivariate distributions.
  \emph{Journal of the American Statistical Association} \bold{83}, 834--841.

  Nelsen, R. B. (2006)
  \emph{An introduction to Copulas}. Springer, New York.
}
\seealso{
  the \code{\linkS4class{copula}} and \code{\linkS4class{acopula}}
  classes, the acopula families, \code{\link{acopula-families}}.
  Constructor functions such as
  \code{\link{ellipCopula}}, \code{\link{archmCopula}},
  \code{\link{fgmCopula}}.
}
\examples{
norm.cop <- normalCopula(0.5)
norm.cop
x <- rcopula(norm.cop, 100)
plot(x)
dcopula(norm.cop, x)
pcopula(norm.cop, x)
persp(norm.cop, dcopula)
contour(norm.cop, pcopula)

## a 3-dimensional normal copula
u <- rcopula(normalCopula(0.5, dim = 3), 1000)
if(require("scatterplot3d"))
  scatterplot3d(u)

## a 3-dimensional clayton copula
cl3 <- claytonCopula(2, dim = 3)
v <- rcopula(cl3, 1000)
pairs(v)
if(require("scatterplot3d"))
  scatterplot3d(v)

## Compare with the "nacopula" version :
fu1 <- dcopula(cl3, v)
fu2 <- copClayton@dacopula(v, theta = 2)
Fu1 <- pcopula(cl3, v)
Fu2 <- pnacopula(onacopula("Clayton", C(2.0, 1:3)), v)
## The density and cumulative values are the same:
stopifnot(all.equal(fu1, fu2, tol= 1e-14),
          all.equal(Fu1, Fu2, tol= 1e-15))
}
\keyword{distribution}
\keyword{multivariate}
