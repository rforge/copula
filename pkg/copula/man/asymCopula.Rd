\name{asymCopula}
\title{Asymmetric Copula Constructors}
\alias{asymBivCopula}
\alias{asymExplicitCopula}
\description{
  Creates an asymmetric copula, an \R object either of class
  \code{"\linkS4class{asymBivCopula}"} or class
  \code{"\linkS4class{asymExplicitCopula}"}.
  
  In the bivariate case, constructs
  an asymmetric bivariate copula from two initial bivariate copulas
  \eqn{C} and \eqn{D} using the general principle described in Liebscher
  (2008). The resulting copula is of the form:
  \deqn{C(u_1^{1 - a_1}, u_2^{1 - a_2}) D(u_1^{a_1}, u_2^{a_2})}{
    C(u_1^{1 - a_1}, u_2^{1 - a_2}) D(u_1^{a_1}, u_2^{a_2})}
  where \eqn{a_1} and \eqn{a_2} are \emph{shape parameters} in [0,1].

  The construction known as \emph{Khoudraji's device} (Khoudraji 1995)
  is a special case of the aforementioned one. See also Genest et
  al. (1998).
}
\usage{
asymBivCopula     (copula1, copula2, shapes)
asymExplicitCopula(copula1, copula2, shapes)
}
\arguments{
  \item{copula1, copula2}{each an \code{\linkS4class{copula}}, of the
    same dimension \eqn{d}, bivariate (\eqn{d = 2}) in the case of
    \code{asymBivCopula()}.}
  \item{shapes}{\code{\link{numeric}} vector of length \eqn{d},
    with values
    in \eqn{[0,1]}.}
}
%% \details{
%% }
%\author{Jun Yan and Ivan Kojadinovic}
\value{
  a new \code{"\linkS4class{asymCopula}"} object, an
 \code{"\linkS4class{asymBivCopula}"} or
 \code{"\linkS4class{asymExplicitCopula}"}, respectively.
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
## An asymetric Clayton copula object
ac2 <- asymBivCopula(copula1 = indepCopula(),
                     copula2 = claytonCopula(6),
                     shapes = c(0.4, 0.95))
ac2 ## FIXME: too long 'fullname'
contour(ac2, dCopula, nlevels = 20, main = "dCopula(<asymBivCopula>)")

ac3 <- asymExplicitCopula(copula1 = claytonCopula(2, dim = 3),
                          copula2 = claytonCopula(6, dim = 3),
                          shapes = c(0.2, 0.6, 0.9))
ac3 ## FIXME: too short 'fullname': does not mention 'claytonCopula'
U3 <- rCopula(200, ac3)
pairs2(U3, pch = par("pch"))

## An asymetric Clayton copula object
ac <- asymBivCopula(copula1 = indepCopula(),
                    copula2 = claytonCopula(6),
                    shapes = c(0.4, 0.95))

## Pseudo-random observations from the copula
n <- 1000
u <- rCopula(n, ac)
plot(u)

\dontrun{
## Fitting can be difficult
## Starting values for the parameters are required 
fitCopula(asymBivCopula(copula2 = claytonCopula()),
          start = c(1.1, 0.5, 0.5), data = pobs(u),
          optim.method="Nelder-Mead")

## Goodness-of-fit testing
n <- 100
u <- rCopula(n, ac)

## Warning: Whether the goodness-of-fit tests hold their level
## needs to be tested
gofCopula(asymBivCopula(copula2 = claytonCopula()), pobs(u),
                        start = c(1.1, 0.5, 0.5),
                        optim.method="Nelder-Mead")

gofCopula(asymBivCopula(copula2 = claytonCopula()), pobs(u),
          start = c(1.1, 0.5, 0.5), optim.method="Nelder-Mead",
          sim = "mult")

u <- rCopula(n, gumbelCopula(4))

gofCopula(asymBivCopula(copula2 = claytonCopula()), pobs(u),
          start = c(1.1, 0.5, 0.5), optim.method="Nelder-Mead")

gofCopula(asymBivCopula(copula2 = claytonCopula()), pobs(u),
          start = c(1.1, 0.5, 0.5), optim.method="Nelder-Mead",
          sim = "mult")
} 

## A "nested" asymetric copula object
agac <- asymBivCopula(copula1 = gumbelCopula(3),
                      copula2 = ac,
                      shapes = c(0.4, 0.95))

}
\keyword{distribution}
\keyword{multivariate}
