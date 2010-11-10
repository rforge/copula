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
  \item{n}{integer specifying the sample size to be used for the random
    variates \eqn{V_0}{V0} and \eqn{V_{01}{V01}}.}
  \item{family}{the Archimedean family (class
    \code{"\linkS4class{acopula}"}) for which \eqn{V_0}{V0} and
    \eqn{V_{01}{V01} are sampled}.}
  \item{taus}{numeric vector of Kendall's taus. It is converted to a vector of
    copula parameters \eqn{theta}, which then serve as
    \eqn{theta_0}{theta0} and \eqn{theta_1}{theta1} for a nested
    Archimedean copula of the specified \code{family}.  First, for each
    \eqn{theta_0}{theta0}, \code{n} random variates \eqn{V_0}{V0} are
    generated.  Then, given the particular \eqn{theta_0}{theta0} and the
    realizations \eqn{V_0}{V0}, \code{n} random variates
    \eqn{V_{01}}{V01} are generated for each \eqn{theta_1}{theta1}
    fulfilling the sufficient nesting condition, see \code{paraConstr}
    in \code{\linkS4class{acopula}}.}
  \item{digits}{number of digits for the output and matrix \code{\link{dimnames}}.}
  \item{verbose}{logical indicating if timing output should generated
    while the random variates are generated (defaults to \code{FALSE}).}
}
\value{
  An \eqn{k \times k}{k x k} matrix of timing in milliseconds
  (\code{1000*\link{system.time}(.)[1]}) where \eqn{k} is \code{length(taus)}.
  The first column contains the times for generation of the
  \eqn{V_0}{V0}s. In the further columns, row \eqn{i} (for
  \eqn{{theta_0}_i}{theta0[i]}) contains the times for the
  \eqn{V_{01}}{V01}s for a particular \eqn{theta_0}{theta0} and all the
  admissible \eqn{theta_1}{theta1}s.
}
\author{Marius Hofert}
\seealso{
  The class \code{\linkS4class{acopula}} and our predefined \code{"acopula"}
  family objects in \code{\link{acopula-families}}.
}
\examples{
timing(10000,"Gumbel", taus = c(0.05,(1:9)/10,0.95))
}
\keyword{timing}
