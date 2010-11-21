\name{opower}
\alias{opower}
\title{Outer Power Transformation of Archimedean Copulas}
\usage{
opower(copbase, thetabase)
}
\description{Apply the outer power transformation to a given Archimedean copula.}
\arguments{ 
  \item{copbase}{a "base" copula, i.e., a copula of class "acopula". This copula must be one of the predefined families.}
  \item{thetabase}{the univariate parameter "theta" for the generator of the base copula \code{copbase}. Hence, the copula which is transformed is fixed, i.e., does not depend on a parameter. }
}
\value{a new "acopula" object, namely the outer power copula based on the provided copula family \code{copbase} with fixed parameter \code{thetabase}. The transform introduces a parameter \code{theta}, so one obtains a parametric Archimedean family object as return value.}
\author{Marius Hofert}
\seealso{
  The class \code{\linkS4class{acopula}} and our predefined "acopula"
  family objects in \code{\link{acopula-families}}.
}
\references{
	Hofert, M. (2010),
	  \emph{Sampling Nested Archimedean Copulas with Applications to CDO Pricing},
	  Suedwestdeutscher Verlag fuer Hochschulschriften AG & Co. KG.

	Hofert, M. (2011),
	  Efficiently sampling nested Archimedean copulas,
	  \emph{Computational Statistics & Data Analysis}, \bold{55}, 57--70.
}
\examples{
## Construct a bivariate Clayton copula with parameter thetabase
thetabase <- .5
C2 <- onacopula("C", C(thetabase, 1:2))

## Compute the corresponding lower and upper tail-dependence coefficients
copClayton@lambdaL(thetabase)
copClayton@lambdaU(thetabase) # => no upper tail dependence

## Now construct a bivariate outer power Clayton copula
opClayton <- opower(copClayton, thetabase)

## Evaluate the tail-dependence coefficients for this copula for the 
## given theta
theta <- 1.5
opClayton@lambdaL(theta)
opClayton@lambdaU(theta) # => upper tail dependence

## Based on opClayton, sample n vectors of a three-dimensional nested 
## outer power Clayton copula with parameters chosen such that the 
## Kendall's tau of the respective bivariate margins are 0.4 and 0.6.
n <- 500
theta0 <- opClayton@tauInv(.4)
theta1 <- opClayton@tauInv(.6) 
## (1) draw V0 and V01
V0  <- opClayton@ V0(n, theta0)
V01 <- opClayton@V01(V0, theta0, theta1)
## (2) build U
U <- cbind(
opClayton@psi(rexp(n)/V0,  theta0),
opClayton@psi(rexp(n)/V01, theta1),
opClayton@psi(rexp(n)/V01, theta1))

## Plot the generated vectors of random variates of the nested outer 
## power Clayton copula.
splom2(U)
}
\keyword{outer power transformation}
