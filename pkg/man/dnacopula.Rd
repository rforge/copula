\name{dnacopula}
\alias{dnacopula}
\title{Copula Density Evaluation}
\description{
  Evaluates the density of an Archimedean copula.
}
\usage{
dnacopula(x, u, MC, log = FALSE)
}
\arguments{
  \item{x}{an object of class \code{"\linkS4class{outer_nacopula}"}.}
  \item{u}{argument of the copula x. Note that u can be a matrix in which case
			  the density is computed for each row of the matrix and the vector of
			  values is returned.}
  \item{MC}{(optional) integer or \code{NULL}. If specified (not NULL)
    Monte Carlo is used for evaluation of the density with sample size
    equal to \code{MC}.}
  \item{log}{if TRUE the log-density is evaluated (defaults to FALSE).}
}
\details{
	The density of an Archimedean copula \eqn{C} (if it exists) with generator 
	\eqn{\psi}{psi} at \eqn{\bm{u}}{u} is given by \deqn{C(\bm{u})=\frac{\psi^{(d)}(\psi^{-1}(u_1)+\dots+\psi^{-1}(u_d))}{\prod_{j=1}^d\psi^\prime(\psi^{-1}(u_j))},\ \bm{u}\in(0,1)^d.}{C(u)=psi^{(d)}(psi^{-1}(u_1)+...+psi^{-1}(u_d))/(psi^prime(psi^{-1}(u_1))*...*psi^prime(psi^{-1}(u_d))), u in (0,1)^d.}
}
\value{
  	A \code{\link{numeric}} vector containing the values of the density of the
	Archimedean copula at \code{u}.
}
\references{
	Hofert, M. and  McNeil, A. J. (2011), to be submitted.	
}
\seealso{
 	For more details about the derivatives of an Archimedean generator, see, e.g., \code{psiDabs} in \code{\link{acopula-class}}.
}
\examples{
## Construct a twenty-dimensional Gumbel copula with parameter chosen 
## such that Kendall's tau of the bivariate margins is 0.25.
theta <- copJoe@tauInv(.25)
C20 <- onacopula("J", C(theta, 1:20))

## Evaluate the copula density at the point u = (0.5,...,0.5)
u <- rep(0.5, 20)
dnacopula(C20, u)

## the same with Monte Carlo based on 10000 simulated "frailties"
dnacopula(C20, u, MC = 10000)

## Evaluate the exact log-density at several points
u <- matrix(runif(100), ncol=20)
dnacopula(C20, u, log = TRUE)	
}
\keyword{distribution}

