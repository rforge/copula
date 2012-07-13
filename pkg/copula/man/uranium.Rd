\name{uranium}
\Rdversion{1.1}
\alias{uranium}
\docType{data}
\title{Uranium Exploration Dataset of Cook & Johnson (1986)}
\description{
  These data consist of log concentrations of 7 chemical elements in 655
  water samples collected near Grand Junction, CO (from the Montrose
  quad-rangle of Western Colorado).  Concentrations were
  measured for the following elements: uranium (U), lithium (Li), cobalt
  (Co), potassium (K), cesium (Cs), scandium (Sc), and titanium (Ti).
}
\usage{data(uranium)}
\format{
  A data frame with 655 observations of the following 7 variables:
  \describe{
    \item{\code{U}}{(numeric) log concentration of uranium.}
    \item{\code{Li}}{(numeric) log concentration of lithium.}
    \item{\code{Co}}{(numeric) log concentration of colbalt.}
    \item{\code{K}}{(numeric) log concentration of potassium.}
    \item{\code{Cs}}{(numeric) log concentration of cesium.}
    \item{\code{Sc}}{(numeric) log concentration of scandum.}
    \item{\code{Ti}}{(numeric) log concentration of titanium.}
  }
}
\references{
  Cook, R. D.  and  Johnson, M. E. (1986)
  Generalized BurrParetologistic distributions with applications to a
  uranium exploration data set.
  \emph{Technometrics} \bold{28}, 123--131.
}
\examples{
data(uranium)
}
\keyword{datasets}
