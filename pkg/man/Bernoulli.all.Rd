\name{Bernoulli.all}
\alias{Bernoulli.all}
\title{Bernoulli Numbers up to Order n}
\description{
  Generate all Bernoulli numbers up to the n-th,
  using the Akiyama-Tanigawa algorithm.
}
\usage{
Bernoulli.all(n, verbose=getOption("verbose"))
}
\arguments{
  \item{n}{positive integer, indicating the index of the largest (and
    last) of the Bernoulli numbers needed.}
  \item{precBits}{\code{NULL} or a positive integer indicating the
    precision of the initial numbrs in bits, using \pkg{"Rmpfr"}'s
    package multiprecision arithmetic.}
  \item{verbose}{logical indicating if the intermediate results of the
    algorithm should be printed.}
}
\value{numeric vector of length n, containing B(n)}
\author{Martin Maechler, 25 Jun 2011 (night train Vienna - Zurich).
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\references{
  Kaneko, Masanobu (2000)
  The Akiyama-Tanigawa algorithm for Bernoulli numbers;
  Journal of Integer Sequences \bold{3}, article 00.2.9
}
\seealso{
  \code{\link{Stirling}}, etc.
}
\examples{
## The example for the paper
MASS::fractions(Bernoulli.all(8, verbose=TRUE))

B10 <- Bernoulli.all(10)
MASS::fractions(B10)

system.time(B50  <- Bernoulli.all(50))# still "no time"
system.time(B100 <- Bernoulli.all(100))# still less than a milli second

## BUT -- the algorithm is *really* not accurate enough, using double prec:
system.time(B100.250 <- as.numeric(Bernoulli.all(100, prec = 250)))
## 0.75 sec [Core i5 (2010)]
m <- cbind(Bn = B100.250, "log10(rel.Err)" =
           round(log1p( - B100/B100.250)/log(10), 2))
rownames(m) <- paste("n=",0:100, sep="")
m[1:5,]
m[2*(1:15) -1,] ## for n=10: still 8 correct digits

system.time(B100.1k <- as.numeric(Bernoulli.all(100, prec = 1024)))
## The first 34 are "the same", but after [41],
## even 250 precBits were *not* sufficient:
round(log10(abs(1 - B100.250/B100.1k))[seq(1,99,by=2)], 2)
}
\keyword{arithmetic}
