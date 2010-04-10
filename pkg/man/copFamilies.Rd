\name{copFamilies}
\title{Specific "ACopula" Family Objects}
\concept{copula Families}
\concept{Archimedean copula family}
\alias{ACopula-families}
\alias{copAMH}
\alias{copClayton}
\alias{copFrank}
\alias{copGumbel}
\alias{copJoe}
\docType{data}
\description{
  Specific \code{"\linkS4class{ACopula}"} family objects are part of
  package \pkg{nacopula}.
}
\usage{
copAMH
copClayton
copFrank
copGumbel
copJoe
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

  \code{copClayton}: ....

  \code{copFrank}: ....

  \code{copGumbel}: ....

  \code{copJoe}: ....

}
\references{
Marius Hofert
}
\examples{
%% FIXME : use something more interesting / reasonable
%% ----   keep at least one example per family
copAMH
copClayton
copFrank
copGumbel
copJoe
}
\keyword{datasets}
