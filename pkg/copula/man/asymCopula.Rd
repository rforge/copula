\name{khoudrajiCopula}
\title{Construction of copulas using Khoudraji's device}
\alias{khoudrajiCopula}
\description{

  Creates an object representing a copula constructed using
  \emph{Khoudraji's device} (Khoudraji, 1995). The resulting \R object
  is either of class \code{"\linkS4class{khoudrajiBivCopula}"} or class
  \code{"\linkS4class{khoudrajiExplicitCopula}"}.

  In the bivariate case, given two copulas \eqn{C_1} and \eqn{C_2},
  Khoudraji's device consists of defining a copula whose c.d.f. is given
  by:

  \deqn{C_1(u_1^{1 - a_1}, u_2^{1 - a_2}) C_2(u_1^{a_1}, u_2^{a_2})}{
    C_1(u_1^{1 - a_1}, u_2^{1 - a_2}) C_2(u_1^{a_1}, u_2^{a_2})}

  where \eqn{a_1} and \eqn{a_2} are \emph{shape parameters} in [0,1].

  The construction principle (see also Genest et al. 1998) is a special
  case of that considered in Liebscher (2008).

}
\usage{
khoudrajiCopula(copula1 = indepCopula(), copula2 = indepCopula(),
                shapes = c(1,1))
}
\arguments{
  \item{copula1, copula2}{each an \code{\linkS4class{copula}}, of the
    same dimension \eqn{d}.}
  \item{shapes}{\code{\link{numeric}} vector of length \eqn{d},
    with values in \eqn{[0,1]}.}
}
%% \details{
%% }
%\author{Jun Yan and Ivan Kojadinovic}
\value{
  A new object of class \code{"\linkS4class{khoudrajiBivCopula}"} in
  dimension two or of class
  \code{"\linkS4class{khoudrajiExplicitCopula}"} when \eqn{d > 2}.
}
\references{
  Genest, C., Ghoudi, K., and Rivest, L.-P. (1998), Discussion of
  "Understanding relationships using copulas", by Frees, E., and
  Valdez, E., \emph{North American Actuarial Journal} \bold{3},
  143--149.

  Khoudraji, A. (1995), Contributions à l'\enc{é}{e}tude des copules et
  \enc{à}{a}la mod\enc{é}{e}lisation des valeurs extr\enc{ê}{e}mes
  bivari\enc{é}{e}es, \emph{PhD thesis, Université Laval}, Québec,
  Canada.

  Liebscher, E. (2008), Construction of asymmetric multivariate copulas,
  \emph{Journal of Multivariate Analysis} \bold{99}, 2234--2250.
}
%% \seealso{
%% }
\examples{
## A bivariate Khoudraji-Clayton copula
kc <- khoudrajiCopula(copula2 = claytonCopula(6),
                      shapes = c(0.4, 0.95))
kc ## FIXME: too long 'fullname'
contour(kc, dCopula, nlevels = 20, main = "dCopula(<khoudrajiBivCopula>)")

## A Khoudraji-Clayton copula with second shape parameter fixed
kcf <- khoudrajiCopula(copula2 = claytonCopula(6),
                       shapes = fixParam(c(0.4, 0.95), c(FALSE, TRUE)))

## A "nested" Khoudraji bivariate copula
kgkcf <- khoudrajiCopula(copula1 = gumbelCopula(3),
                         copula2 = kcf,
                         shapes = c(0.7, 0.25))
contour(kgkcf, dCopula, nlevels = 20, main = "dCopula(<khoudrajiBivCopula>)")

n <- 300
u <- rCopula(n, kc)
plot(u)

\dontrun{
## Fitting can be difficult
## Starting values are required for all parameters
fitCopula(khoudrajiCopula(copula2 = claytonCopula()),
          start = c(1.1, 0.5, 0.5), data = pobs(u),
          optim.method = "Nelder-Mead")

## Second shape parameter fixed to 1
fitCopula(kcf,
          start = c(1.1, 0.5), data = pobs(u),
          optim.method = "Nelder-Mead")

## With a different optimization method
fitCopula(kcf,
          start = c(1.1, 0.5), data = pobs(u),
          optim.method = "BFGS")

## GOF example
gofCopula(kcf, x = u, start = c(1.1, 0.5), optim.method = "BFGS")
gofCopula(kcf, x = u, start = c(1.1, 0.5), optim.method = "BFGS",
          sim = "mult")
## Warning: Whether the goodness-of-fit tests hold their level
## needs to be tested

## Another example
u <- rCopula(n, gumbelCopula(4))
gofCopula(kcf, x = u, start = c(1.1, 0.5), optim.method = "BFGS")
gofCopula(kcf, x = u, start = c(1.1, 0.5), optim.method = "BFGS",
          sim = "mult")
}
}
\keyword{distribution}
\keyword{multivariate}
