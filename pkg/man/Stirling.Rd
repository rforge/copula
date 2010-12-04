\name{Stirling}
\title{Stirling Numbers of First and Second Kind}
\alias{Stirling1}
\alias{Stirling2}
\alias{Stirling1.all}
\alias{Stirling2.all}
\description{
  Compute Stirling numbers of the first and second kind,
  possibly for all \eqn{k} \dQuote{at once}.
}
\usage{
Stirling1(n, k)
Stirling2(n, k, method = c("lookup.or.store", "direct"))

Stirling1.all(n)
Stirling2.all(n)
}
\arguments{
  \item{n}{positive integer}
  \item{k}{integer in \code{0:n}.}
  \item{method}{for \code{Stirling2()}, string specifying the method to
    be used.  \code{"direct"} uses the explicit formula.}
}
\details{
  Stirling numbers of the first kind:\cr
  s(n,k) = (-1)^{n-k} times
  the number of permutations of 1,2,\dots,n with exactly k cycles.

  Stirling numbers of the second kind:\cr
  \eqn{S^{(k)}_n}{S(n,k)} is the number of ways of partitioning a set
    of \eqn{n} elements into \eqn{k} non-empty subsets.
}
\value{
  \eqn{s(n,k)} or \eqn{S(n,k) = S^{(k)}_n}{S(n,k)}, respectively.

  \code{Stirling1.all(n)} is the same as \code{sapply(1:n, Stirling1, n=n)},
  and\cr
  \code{Stirling2.all(n)} is the same as \code{sapply(1:n, Stirling2, n=n)},
  but more efficiently.
}
\author{Martin Maechler ("direct": May 1992)}
\references{
  Abramowitz and Stegun
  24,1,4 (p. 824-5 ; Table 24.4, p.835);
  Closed Form : p.824 "C."

  NIST Digital Library of Mathematical Functions,
  26.8: \url{http://dlmf.nist.gov/26.8}
}
% \seealso{
% }
\examples{
Stirling1(7,2)
Stirling2(7,3)

Stirling1.all(9)
Stirling2.all(9)
}
\keyword{arithmetic}
