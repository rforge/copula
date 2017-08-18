\name{margCopula}
\title{Marginal copula of a copula with specified margins}
\alias{margCopula}
\description{
  The marginal copula of a copula is 
}
\usage{
martCopula(copula, keep)
}
\arguments{
  \item{copula}{object of class \code{"\linkS4class{copula}"}}
  \item{keep}{logical vector indicating which margins to keep}
}
\details{
  The marginal copula of a copula is needed in practical data analysis when
  one or more of the components of some multivariate observations is
  missing. For normal/t/Archimedean copulas, the marginal copulas can be
  easily obtained. For a general copula, this may not be an easy problem.
  The current implementation only supports normal/t/Archimedean copulas.
}
\value{
  The marginal copula of the specified margin.
}
%\author{Jun Yan}
\examples{
tc <- tCopula(8:2 / 10, dim = 8, dispstr = "toep")
margCopula(tc, c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE))

nc <- normalCopula(.8, dim = 8, dispstr = "ar1")
margCopula(nc, c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE))

gc <- gumbelCopula(2, dim = 8)
margCopula(gc, c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
}
\keyword{models}
\keyword{multivariate}
