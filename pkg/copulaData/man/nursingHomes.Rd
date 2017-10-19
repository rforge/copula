\name{nursingHomes}
\alias{nursingHomes}
\title{Wisconsin Nursing Homes Utilization Data}
\description{
  Data set containing the occupancy rate (utilization) and
  covariates of 377 nursing homes in Wisconsin between 1995 and 2001.
}
\usage{
data("nursingHomes")
}
\format{
  \code{\link{data.frame}} containing the nursing home identification
  number (\code{ID}), the occupancy rate (\code{Rate}; see Sun et
  al. (2008, Equation (7)) for how this is computed), the logarithmic
  number of beds of the nursing home (\code{LnNumBed}), the logarithmic
  net square foot of the nursing home (\code{LnSqrFoot}), the cost
  report year (\code{CRYear}), indicators whether the nursing home runs
  on a for-profit basis (\code{Pro}), whether it is tax exempt
  (\code{TaxExempt}), whether it has self funding of insurance
  (\code{SelfIns}), whether it is accredited as Medicare Certified
  (\code{MCert}) and whether it is located in an urban environment
  (\code{Urban}); see Sun et al. (2008, Table 2).
}
\source{
  Edward W. Frees,
  Wisconsin Department of Health and Family Services (by now named
  \dQuote{Wisconsin Department of Health Services})
}
\references{
  Sun, J., Frees, E. W. and Rosenberg, M. A. (2008)
  Heavy-tailed longitudinal data modeling using copulas.
  \emph{Insurance: Mathematics and Economics}
  \bold{42}, 817--830.
}
\examples{
data("nursingHomes")
str(nursingHomes)
}
\keyword{datasets}