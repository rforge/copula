\name{mixCopula}
\alias{mixCopula}
\title{Create Mixture of Copulas}
\description{
  A mixture of \eqn{m} copulas of dimension \eqn{d} with weights
  \eqn{w_j}, \eqn{j=1,2,\ldots,m} is itself a \eqn{d}-dimensional copula,
  with cumulative distribution function
  \deqn{C(x) = \sum_{j=1}^m w_j C_j(x),}{C(x) = sum(j=1..m; w[j] C[j]),}
  and (probability) density function
  \deqn{c(x) = \sum_{j=1}^m w_j c_j(x),}{c(x) = sum(j=1..m; w[j] c[j]),}
  where \eqn{C_j}{C[j]} are the CDFs and \eqn{c_j}{c[j]} are the
  densities of the \eqn{m} component copulas,
  \eqn{j=1,2,\ldots,m}.
}
\usage{
mixCopula(coplist, w = NULL)
}
\arguments{
  \item{coplist}{a \code{\link{list}} of length \eqn{m (\ge 1)}{m (>=
      1)} copulas (each inheriting from \code{\linkS4class{parCopula}}), all
    of the same dimension.}
  \item{w}{numeric vector of length \eqn{m} of non-negative mixture
    weights, or \code{NULL}, which means equal weights.}
}
\value{
  an object of \code{\link{class}} \code{\linkS4class{mixCopula}}
}
\details{
  It easy to see that the tail dependencies \code{\link{lambda}()} and
  Spearman's rank correlation \code{\link{rho}()} can be computed as
  mixture of the individual measures.% but Kendall's tau cannot !
}
% \author{ Martin Maechler }
\seealso{
 \code{\link{khoudrajiCopula}}, \code{\link{rotCopula}} also create new
 copula models from existing ones.
}
\examples{%% NB: Now also have ../tests/mixCop-tst.R
mC <- mixCopula(list(gumbelCopula(2.5, dim=3),
                     claytonCopula(pi, dim=3),
                     tCopula(0.7, dim=3)),
                c(2,2,4)/8)
mC
stopifnot(dim(mC) == 3)

set.seed(17)
uM <- rCopula(600, mC)
splom2(uM, main = "mixCopula( (gumbel, clayton, t-Cop) )")
d.uM <- dCopula(uM, mC)
p.uM <- pCopula(uM, mC)

## mix a Gumbel with a rotated Gumbel (with equal weights 1/2):
mGG <- mixCopula(list(gumbelCopula(2), rotCopula(gumbelCopula(1.5))))
rho(mGG)   # 0.57886
lambda(mGG)# both lower and upper tail dependency

loglikCopula(c(2.5, pi, rho.1=0.7, df = 4, w = c(2,2,4)/8),
                    u = uM, copula = mC)
ll.df <- Vectorize(function(df, rho)
                   loglikCopula(c(2.5, pi, rho.1=rho, df=df, w = c(2,2,4)/8),
                                uM, mC))
(df. <- 1/rev(seq(1/8, 1/2, length=21)))# in [2, 8] equidistant in 1/. scale

ll. <- ll.df(df., rho = 0.7)
plot(df., ll., type = "b", ylim = extendrange(ll., f = 0.10),
     main = "loglikCopula((.,.,., df, .,.,.), u, <mixCopula>)")

if(!exists("Xtras")) Xtras <- copula:::doExtras() ; cat("Xtras: ", Xtras,"\n")
if (Xtras) {
%%% FIXME: use a grid of 'Rhos' and then use  persp() or wireframe2() 
  Rhos <- c(0.7, 0.75, 0.65, 0.72, 0.78, 0.77)
  Cols <- c(palette()[c(1,2,4)], "tomato", "forest green", "dark green")
  for(k in 2:length(Rhos))
    lines(df., ll.df(df., rho = Rhos[k]), type = "b", col=Cols[k])
  legend("bottom", pch = 1, lwd = 1, col = Cols, bty = "n",
       as.expression(lapply(Rhos, function(.) substitute(rho == R., list(R. = .)))))
}

## TODO: fitCopula() example with "fixed weights"


}% {examples}
\keyword{distribution}
\keyword{multivariate}
