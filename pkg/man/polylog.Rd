\name{polylog}
\alias{polylog}
\title{Polylogarithm Li_s(z)}
\description{
  Compute the polylogarithm function \eqn{Li_s(z)},
  initially defined as the power series,
  \deqn{\mathop{Li}_s(z) = \sum_{k=1}^\infty {z^k \over k^s},}{%
    Li_s(z) = sum(k=1..Inf; z^k / k^s),}
  for \eqn{|z| < 1}, and then more generally (by analytic continuation) as
  \deqn{\mathop{Li}_1(z) = -\mathop{ln}(1-z),}{Li_1(z) = -ln(1-z),}
    and
  \deqn{\mathop{Li}_{s+1}(z) = \int_0^z \frac{\mathop{Li}_s(t)}{t}\,dt.}{%
    Li_{s+1}(z) = Int[0..z] (Li_s(t) / t) dt.}

  Currently only the case of negative integer \eqn{s} is well supported,
  as that is used for some of the Archimedean copula densities.
}
\usage{
polylog(z, s,
        method = c("sum", "negI-s-Stirling",
                   "negI-s-Eulerian", "negI-s-asymp-w"),
        logarithm = FALSE, is.log.z = FALSE, is.logmlog = FALSE,
        asymp.w.order = 0, n.sum)
}
\arguments{
  \item{z}{numeric or complex vector}
  \item{s}{complex number; current implementation is aimed at
    \eqn{s \in \{0,-1,\dots}{s in 0,-1,...}.}
  \item{method}{a string specifying the algorithm to be used.}
  \item{logarithm}{logical specifying to return log(Li.(.)) instead of Li.(.).}
  \item{is.log.z}{logical; if TRUE, the specified \code{z} argument is really  w = log(z);
    i.e., compute  \eqn{Li_s(exp(w))}, and we typically have \eqn{w <
      0}, or equivalently, \eqn{z < 1}.}
  \item{is.logmlog}{logical; if TRUE, the specified \code{z} argument is
    really .... FIXME...} % log(-log(z)) w = log(z);
  \item{asymp.w.order}{only default is implemented currently...}
  \item{n.sum}{for \code{method="sum"} only: the number of terms used.}
}
\details{
  Almost entirely taken from
  \url{http://en.wikipedia.org/wiki/Polylogarithm}:

For integer values of the polylogarithm order, the following
explicit expressions are obtained by repeated application of
\eqn{z \cdot \frac{\partial}{\partial z}}{z * d/dz} to \eqn{Li_1(z)}:
%---
\deqn{
  {Li}_{1}(z) = -\ln(1-z) \\
  {Li}_{0}(z) = {z \over 1-z} \\
  {Li}_{-1}(z) = {z \over (1-z)^2} \\
  {Li}_{-2}(z) = {z \,(1+z) \over (1-z)^3} \\
  {Li}_{-3}(z) = {z \,(1+4z+z^2) \over (1-z)^4}
}{%
  {Li}_{1}(z) = -ln(1-z),
  {Li}_{0}(z) =  z / (1-z),
  {Li}_{-1}(z) = z / (1-z)^2,
  {Li}_{-2}(z) = z (1+z) / (1-z)^3,
  {Li}_{-3}(z) = z (1+4z+z^2) / (1-z)^4.
}

Accordingly, the polylogarithm reduces to a ratio of polynomials in
z, and is therefore a rational function of z, for all nonpositive
integer orders. The general case may be expressed as a finite sum:
%---

\deqn{\mathop{Li}_{-n}(z) =
  \left(z \,{\partial \over \partial z} \right)^n \frac{z}{1-z} =
  = \sum_{k=0}^n k! \,S(n+1,k+1) \left({z \over {1-z}} \right)^{k+1}
  \ \ (n=0,1,2,\ldots),}{%
  Li[-n](z) = ( z d/dz )^n  z/(1-z) =
            = sum(k=0..n ; k! S(n+1,k+1) (z /(1-z))^(k+1)),   (n=0,1,2,...),}
where \eqn{S(n,k)} are the Stirling numbers of the second kind.

Equivalent formulae applicable to negative integer orders are
(Wood 1992, § 6) ...
\deqn{\mathop{Li}_{-n}(z) = {1 \over (1-z)^{n+1}} \sum_{k=0}^{n-1} \left\langle {n \atop k} \right\rangle z^{n-k} \qquad (n=1,2,3,\ldots) ~, }{%
  {Li}_{-n}(z) = \frac{1}{(1-z)^{n+1}} sum_{k=0}^{n-1} < n \ k >  z^{n-k}
  = \frac{z \sum_{k=0}^{n-1} < n \ k >  z^k} {(1-z)^{n+1}},}
where  \eqn{\left\langle {n \atop k} \right\rangle}{< n \ k >}  are the
Eulerian numbers, see also \code{\link{Eulerian}}.

% All roots of Li_{-n}(z) are distinct and real; they include z = 0.

% Duplication formula:  2^{1-s} Li_s(z^2) = Li_s(z) + Li_s(-z).
}
\value{
  numeric/complex vector as \code{z}.
}
\author{Martin Maechler}
\references{
  Wikipedia (2011) \emph{Polylogarithm},
  \url{http://en.wikipedia.org/wiki/Polylogarithm}.

  Wood, D.C. (June 1992).
  The Computation of Polylogarithms. Technical Report 15-92.
  Canterbury, UK: University of Kent Computing Laboratory.
  \url{http://www.cs.kent.ac.uk/pubs/1992/110}.% Retrieved 2005-11-01.

  Apostol, T.M. (2010), \emph{"Polylogarithm"}, in the
  NIST Handbook of Mathematical Functions, \url{http://dlmf.nist.gov/25.12}

  Lewin, L. (1981).
  \emph{Polylogarithms and Associated Functions}.
  New York: North-Holland. ISBN 0-444-00550-1.
}
\seealso{
  is used in MLE for some A.copulas, see \code{\link{emle}}.
}
\examples{
polylog(z = 1, s = 2, n.sum = 1e5)
## in the limit, should be equal
pi^2 / 6

p.Li <- function(s.set, from = -2.6, to = 1/4, ylim = c(-1, 0.5),
                 colors = c("orange","brown", palette()), n = 201, ...)
{
    s.set <- sort(s.set, decreasing = TRUE)
    s <- s.set[1] # <_ for auto-ylab
    curve(polylog(x, s, method="negI-s-Stirling"), from, to,
          col=colors[1], ylim=ylim, n=n, ...)
    abline(h=0,v=0, col="gray")
    for(is in seq_along(s.set)[-1])
        curve(polylog(x, s=s.set[is], method="negI-s-Stirling"),
              add=TRUE, col = colors[is], n=n)
    s <- rev(s.set)
    legend("bottomright", paste("s =", s), col=colors[2-s], lty=1, bty="n")
}

## yellow is unbearable (on white):
palette(local({p <- palette(); p[p=="yellow"] <- "goldenrod"; p}))

## Wikipedia page plot (+/-):
p.Li(1:-3, ylim= c(-.8, 0.6), colors = c(2:4,6:7))

## and a bit more:
p.Li(1:-5)

## For the range we need it:
ccol <- c(NA,NA, rep(palette(),10))
p.Li(-1:-20, from=0, to=.99, colors=ccol, ylim = c(0, 10))
## log-y scale:
p.Li(-1:-20, from=0, to=.99, colors=ccol, ylim = c(.01, 1e7), log="y", yaxt="n")
if(require("sfsmisc")) eaxis(2) else axis(2)
}
\keyword{arithmetic}
