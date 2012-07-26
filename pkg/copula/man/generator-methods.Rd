\name{generator}
\docType{methods}
%----------------------Archimedean Copulas------------------
% iPsi
\alias{iPsi}
\alias{iPsi-methods}
\alias{iPsi,claytonCopula-method}
\alias{iPsi,frankCopula-method}
\alias{iPsi,gumbelCopula-method}
\alias{iPsi,amhCopula-method}
\alias{iPsi,joeCopula-method}
% psi
\alias{genInv}%deprecated
\alias{psi}
\alias{psi-methods}
\alias{psi,claytonCopula-method}
\alias{psi,frankCopula-method}
\alias{psi,gumbelCopula-method}
\alias{psi,amhCopula-method}
\alias{psi,joeCopula-method}
% iPsi
\alias{genFun}%deprecated
\alias{iPsi}
\alias{iPsi-methods}
\alias{iPsi,claytonCopula-method}
\alias{iPsi,frankCopula-method}
\alias{iPsi,gumbelCopula-method}
\alias{iPsi,amhCopula-method}
\alias{iPsi,joeCopula-method}
% diPsi
\alias{genFunDer1}%deprecated
\alias{genFunDer2}%deprecated
\alias{diPsi}
\alias{diPsi-methods}
\alias{diPsi,claytonCopula-method}
\alias{diPsi,frankCopula-method}
\alias{diPsi,gumbelCopula-method}
\alias{diPsi,amhCopula-method}
\alias{diPsi,joeCopula-method}
%----------------------Extreme Value Copulas------------------
% A
\alias{Afun}%deprecated
\alias{A}
\alias{A-methods}
\alias{A,galambosCopula-method}
\alias{A,gumbelCopula-method}
\alias{A,huslerReissCopula-method}
\alias{A,tawnCopula-method}
\alias{A,tevCopula-method}
\alias{A,asymCopula-method}
\alias{A,indepCopula-method}
% dAdu
\alias{AfunDer}%deprecated
\alias{dAdu}
\alias{dAdu-methods}
\alias{dAdu,galambosCopula-method}
\alias{dAdu,gumbelCopula-method}
\alias{dAdu,huslerReissCopula-method}
\alias{dAdu,tawnCopula-method}
\alias{dAdu,tevCopula-method}
% main
\title{Generator functions for Archimedean and extreme value copulas}
\description{
  Methods to evaluate the generator function, the inverse generator
  function, and derivatives of the generator function.
}
\usage{
psi(copula, s)
%% TODO?  dPsi(copula, s)
iPsi(copula, u, \dots)
diPsi(copula, u, degree=1, log=FALSE, \dots)

A(copula, w)
dAdu(copula, w)
}
\arguments{
  \item{copula}{an object of class \code{"\linkS4class{copula}"}. }
  \item{u, s, w}{numerical vector at which these functions are to be
    evaluated. }
  \item{\dots}{further arguments for specific families.}
  \item{degree}{the degree of the derivative (defaults to 1).}
  \item{log}{logical indicating if the \code{\link{log}} of the
    \emph{absolute} derivative is desired.  Note that the derivatives
    alternate in sign.}
}
\details{
  \code{psi()} and \code{iPsi()} are, respectively, the generator
  function and its inverse for an Archimedean copula.

  \code{diPsi()} computes (currently only the first two) derivatives of
  the generator function for an Archimedean copula.

  \code{A()} is the generator function of an extreme value copula.

  \code{dAdu()} returns a data.frame containing the 1st and 2nd
  derivative of \code{A()}.
}
\examples{
## List the available methods (and their definitions):
showMethods("A")
showMethods("iPsi", incl=TRUE)
}
\keyword{multivariate}
