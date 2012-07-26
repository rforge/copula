\name{generator}
\docType{methods}
%----------------------Archimedean Copulas------------------
% genFun
\alias{genFun}
\alias{genFun-methods}
\alias{genFun,claytonCopula-method}
\alias{genFun,frankCopula-method}
\alias{genFun,gumbelCopula-method}
\alias{genFun,amhCopula-method}
\alias{genFun,joeCopula-method}
% psi
\alias{genInv}%deprecated
\alias{psi}
\alias{psi-methods}
\alias{psi,claytonCopula-method}
\alias{psi,frankCopula-method}
\alias{psi,gumbelCopula-method}
\alias{psi,amhCopula-method}
\alias{psi,joeCopula-method}
% genFunDer1
\alias{genFunDer1}
\alias{genFunDer1-methods}
\alias{genFunDer1,claytonCopula-method}
\alias{genFunDer1,frankCopula-method}
\alias{genFunDer1,gumbelCopula-method}
\alias{genFunDer1,amhCopula-method}
\alias{genFunDer1,joeCopula-method}
% genFunDer2
\alias{genFunDer2}
\alias{genFunDer2-methods}
\alias{genFunDer2,claytonCopula-method}
\alias{genFunDer2,frankCopula-method}
\alias{genFunDer2,gumbelCopula-method}
\alias{genFunDer2,amhCopula-method}
\alias{genFunDer2,joeCopula-method}
%----------------------Extreme Value Copulas------------------
% Afun
\alias{Afun}
\alias{Afun-methods}
\alias{Afun,galambosCopula-method}
\alias{Afun,gumbelCopula-method}
\alias{Afun,huslerReissCopula-method}
\alias{Afun,tawnCopula-method}
\alias{Afun,tevCopula-method}
\alias{Afun,asymCopula-method}
\alias{Afun,indepCopula-method}
% AfunDer
\alias{AfunDer}
\alias{AfunDer-methods}
\alias{AfunDer,galambosCopula-method}
\alias{AfunDer,gumbelCopula-method}
\alias{AfunDer,huslerReissCopula-method}
\alias{AfunDer,tawnCopula-method}
\alias{AfunDer,tevCopula-method}
% main
\title{Generator functions for Archimedean and extreme value copulas}
\description{
  Methods to evaluate the generator function, the inverse generator
  function, and derivatives of the generator function.
}
\usage{
genFun(copula, u)
psi(copula, s)
iPsi(copula, u, ...)
genFunDer1(copula, u)
genFunDer2(copula, u)
diPsi(copula, u, degree=1, ...)
Afun(copula, w)
AfunDer(copula, w)
}
\arguments{
  \item{copula}{an object of class \code{"\linkS4class{copula}"}. }
  \item{u, s, w}{numerical vector at which these functions are to be
    evaluated. }
}
\details{
  \code{psi()} and \code{genFun()} are, respectively, the generator
  function and its inverse function for an Archimedean copula.

  \code{genFunDer1()} and \code{genFunDer2()} are, respectively, the 1st
  and 2nd derivatives of the generator function for an Archimedean copula.

  \code{Afun()} is the generator function of an extreme value copula.

  \code{AfunDer()} returns a data.frame containing the 1st and 2nd
  derivative of \code{Afun()}.
}
\examples{
## List the available methods (and their definitions):
showMethods("Afun")
showMethods("genFun", incl=TRUE)
}
\keyword{multivariate}
