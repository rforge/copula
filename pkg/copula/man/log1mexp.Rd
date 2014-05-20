\name{log1mexp}
\alias{log1pexp}
\alias{log1mexp}
\title{Compute  f(a) = \eqn{\mathrm{log}}{log}(1 +/- \eqn{\mathrm{exp}}{exp}(-a))
  Numerically Optimally}
\description{
  Compute f(a) = log(1 - exp(-a)), respectively
  g(x) = log(1 + exp(x))  quickly numerically accurately.
}
\usage{
log1mexp(a, cutoff = log(2))
log1pexp(x, c0 = -37, c1 = 18, c2 = 33.3)
}
\arguments{
  \item{a}{numeric vector of positive values}
  \item{x}{numeric vector}
  \item{cutoff}{positive number; \code{log(2)} is \dQuote{optimal},
    %% see below, TODO
    but the exact value is unimportant, and anything in
    \eqn{[0.5, 1]} is fine.}
  \item{c0, c1, c2}{cutoffs for \code{log1pexp}; see below.}
}
\value{
  f(a) == log(1 - exp(-a)) == log1p(-exp(-a)) ==
  log(-expm1(-a))

  or

  g(x) == log(1 + exp(x)) == log1p(exp(x))

  computed accurately and quickly
}
\author{
  Martin Maechler, May 2002; \code{log1pexp()} in 2012.
}
\references{%  ~/R/Pkgs/Rmpfr/vignettes/log1mexp-note.Rnw
  Martin \enc{Mächler}{Maechler} (2012).
  Accurately Computing \eqn{\log(1-\exp(-|a|))};
  \url{http://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf}.
  %  see also <<log1mexp>> in ../inst/doc/Frank-Rmpfr.Rnw
}
\examples{
a <- 2^seq(-58,10, length = 256)
fExpr <- expression(
          log(1 - exp(-a)),
          log(-expm1(-a)),
          log1p(-exp(-a)),
          log1mexp(a))
names(fExpr) <- c("DEF", "expm1", "log1p", "F")
str(fa <- do.call(cbind, as.list(fExpr)))
matplot(a, fa, type = "l", log = "x")# don't see the big difference
legend("topleft", fExpr, col=1:4, lty=1:4, bty="n")

head(fa)# expm1() works here
tail(fa)# log1p() works here

matplot(a, -fa, type = "l", log = "xy")# don't see the big difference
legend("left", fExpr, col=1:4, lty=1:4, bty="n")

curve(log1pexp, -10, 10, asp=1)
abline(0,1, h=0,v=0, lty=3, col="gray")

## Cutoff c1 for log1pexp() -- not often "needed":
curve(log1p(exp(x)) - log1pexp(x), 16, 20, n=2049)
## need for *some* cutoff:
x <- seq(700, 720, by=2)
cbind(x, log1p(exp(x)), log1pexp(x))

## Cutoff c2 for log1pexp():
curve((x+exp(-x)) - x, 20, 40, n=1025)
curve((x+exp(-x)) - x, 33.1, 33.5, n=1025)
}
\keyword{math}
