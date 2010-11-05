\name{timing}
\alias{timing}
\title{Timing for Sampling Nested Archimedean Copulas}
\usage{
timing(n, family, taus, digits = 3, verbose = FALSE)
}
\description{
  This function provides measurements of (user) run times for the frailty
  variables involved in a nested Archimedean copula.
}
\arguments{
	\item{n}{The sample size to be used for the random variates \eqn{V_0}{V0} and
	 	\eqn{V_{01}{V01}}.}
	\item{family}{The Archimedean family for which \eqn{V_0}{V0} and
	 	\eqn{V_{01}{V01} are sampled}.}
	\item{taus}{A vector of Kendall's taus. It is converted to a vector of copula
	 	parameters \eqn{theta}, which then serve as \eqn{theta_0}{theta0} and
	 	\eqn{theta_1}{theta1} for a nested Archimedean copula of the provided
	 	\code{family}. First, for each of the \eqn{theta_0}{theta0}, \code{n}
	 	random variates \eqn{V_0}{V0} are generated. Then, given the particular
	 	\eqn{theta_0}{theta0} and the realizations \eqn{V_0}{V0}, \code{n} random
	 	variates \eqn{V_{01}}{V01} are generated for each \eqn{theta_1}{theta1}
	 	fulfilling the sufficient nesting condition, see \code{paraConstr} in
	 	\code{\linkS4class{acopula}}.}
	\item{digits}{The number of digits for the output (defaults to three).}
	\item{verbose}{If \code{TRUE}, timing output is generated while the random
	 	variates are generated (defaults to \code{FALSE}).}
}
\value{
		An l-by-l matrix where l is the length of \code{taus}. The rows (standing
	   for the choices of \eqn{theta_0}{theta0}) show the timing output for the
	 	\eqn{V_{01}}{V01}s for a particular \eqn{theta_0}{theta0} and all the
	  	admissible \eqn{theta_1}{theta1} (the latter being the columns). The first
	   column shows the timing output for the \eqn{V_0}{V0}s.
}
\author{Marius Hofert}
\seealso{
  The class \code{\linkS4class{acopula}} and our predefined "acopula"
  family objects in \code{\link{acopula-families}}.
}
\examples{
timing(10000,"Gumbel",c(0.05,(1:9)/10,0.95))
}
\keyword{timing}
