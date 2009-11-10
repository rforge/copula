\name{rstable1}
\title{Random Numbers from (Skew) Stable Distribution}
\alias{rstable1}
\alias{rstable}
\description{
  Generate random numbers from the stable distribution,
  using its  \dQuote{S1}- (or \dQuote{1})-parametrization.
}
\usage{
rstable1(n, alpha, beta, gamma = 1, delta = 0)
}
\arguments{
    \item{alpha, beta, gamma, delta}{
        value of the index parameter \code{alpha} with \code{alpha = (0,2]};
        skewness parameter \code{beta}, in the range [-1, 1];
        scale parameter \code{gamma}; and
        shift or location parameter \code{delta}.
        }
    \item{n}{an integer, the number of observations to generate.}
}

\value{
  a \code{\link{numeric}} vector of length \code{n}.
}
\details{
    \bold{Skew Stable Distribution:}

    The function uses the approach of J.P. Nolan for general
    stable distributions.  \code{rstable1} only uses the so called
    \dQuote{S}, \dQuote{S1}, or \dQuote{1} parameterization:
    \cr
    It is the parameterization used by Samorodnitsky and Taqqu in the
    book \emph{Stable Non-Gaussian Random Processes}, and is a slight
    modification of Zolotarev's (A) parameterization.
}
\author{
  Diethelm Wuertz wrote \code{\link[fBasics]{rstable}} for Rmetrics;
  Martin Maechler vectorized it and sped it up.
}
\seealso{\code{\link[fBasics]{rstable}} which allows other
  parametrizations than the \dQuote{S1} one.
}
\references{
  Chambers J.M., Mallows, C.L. and Stuck, B.W. (1976);
  \emph{A Method for Simulating Stable Random Variables},
  J. Amer. Statist. Assoc. 71, 340--344.

  Nolan J.P. (1999);
  \emph{Stable Distributions},
  Preprint, University Washington DC, 30 pages.

% Nolan J.P. (1999);
%     \emph{Numerical Calculation of Stable Densities and Distribution
%     Functions},
%     Preprint, University Washington DC, 16 pages.

% Samoridnitsky G., Taqqu M.S. (1994);
%     \emph{Stable Non-Gaussian Random Processes, Stochastic Models
%     with Infinite Variance},
%     Chapman and Hall, New York, 632 pages.

% Weron, A., Weron R. (1999);
%     \emph{Computer Simulation of Levy alpha-Stable Variables and
%     Processes},
%     Preprint Technical Univeristy of Wroclaw, 13 pages.
}
\examples{
   # Plot rvs Series
   set.seed(1953)
   r <- rstable1(n = 1000, alpha = 1.9, beta = 0.3)
   plot(r, type = "l", main = "stable: alpha=1.9 beta=0.3",
        col = "steelblue")
   grid()
}
\keyword{distribution}

