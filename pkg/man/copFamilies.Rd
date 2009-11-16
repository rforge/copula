\name{copFamilies}
\alias{copAMH}
\alias{copClayton}
\docType{data}
\title{Specific "ACopula" Family Objects}
\description{
  Specific \code{"\linkS4class{ACopula}"} family objects are part of
  package \pkg{nacopula}.
}
\usage{
copAMH
copClayton
}
\format{
  All these are objects of formal class \code{"\linkS4class{ACopula}"}.
  with 13 slots
  ..@ name      : chr "AMH"  (or "Clayton" ....)
  ..@ psi       :function (t, theta)
  ..@ psiInv    :function (t, theta)
  ..@ paraConstr:function (theta)
  ..@ V0        :function (n, theta)
  ..@ tau       :function (theta)
  ..@ tauInv    :function (tau)
  ..@ lTDC      :function (theta)
  ..@ lTDCInv   :function (lambda)
  ..@ uTDC      :function (theta)
  ..@ uTDCInv   :function (lambda)
  ..@ nestConstr:function (theta0, theta1)
  ..@ V01       :function (V0, theta0, theta1)
}
\details{
  \code{copAMH} : ....

  \code{copClayton}:

}
\references{
Marius Hofert
}
\examples{
copAMH
}
\keyword{datasets}
