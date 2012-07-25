## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


joeCopula <- function(param = NA_real_, dim = 2L) {
  new("joeCopula",
      dimension = as.integer(dim),
      parameters = param[1],
      param.names = "param",
      param.lowbnd = 1,
      param.upbnd = Inf,
      message = "Joe copula family; Archimedean copula")
}

setMethod("rcopula", signature("joeCopula"), function(copula, n)
	  racopula(n, Cp=copula, theta = copula@parameters, d = copula@dimension))

setMethod("pcopula", signature("joeCopula"),
	  function (copula, u, ...) pacopula(copJoe, u, theta=copula@parameters))
setMethod("dcopula", signature("joeCopula"),
	  function (copula, u, log = FALSE, ...)
	  copJoe@dacopula(u, theta=copula@parameters, log=log, ...))

setMethod("genFun", signature("joeCopula"),
	  function(copula, u) copJoe@psi(t=u, theta=copula@parameters))
setMethod("genInv", signature("joeCopula"),
	  function(copula, s) copJoe@iPsi(u=s, theta=copula@parameters))
setMethod("genFunDer1", signature("joeCopula"), function(copula, u)
          - copJoe@absdPsi(u, theta=copula@parameters))
setMethod("genFunDer2", signature("joeCopula"), function(copula, u)
          copJoe@absdPsi(u, theta=copula@parameters, degree=2))

setMethod("kendallsTau", signature("joeCopula"),
          function(copula) copJoe@tau(theta=copula@parameters))
setMethod("tailIndex", signature("joeCopula"),
	  function(copula) c(lower=0,
			     upper=copJoe@lambdaU(theta=copula@parameters)))

setMethod("calibKendallsTau", signature("joeCopula"),
	  function(copula, tau) copJoe@tauInv(tau))

## "TODO"
## setMethod("spearmansRho", signature("joeCopula"), ... ? ...)
## setMethod("calibSpearmansRho", signature("joeCopula"), function(copula, rho) ...)

## "TODO"
## setMethod("dRho", signature("joeCopula"), ...)
## setMethod("dTau", signature("joeCopula"), ...)
