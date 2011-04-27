## Copyright (C) 2010 Marius Hofert and Martin Maechler
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

#### List of supported Archimedean copulas

## FIXME: Not "nice" that the names of the copula objects must be be used
##	inside some of their components.
## Possible solution:
##	use *same* environment for  (tau, tauInv, paraConstr, nestConstr)

## NOTA BENE:  Write psi(), tau(), ... functions such that they *vectorize*
## ---------   *and* do work even for (non-numeric) NULL argument
##          	{now checked in "acopula" validityMethod -- see ./AllClass.R }

### ==== Ali-Mikhail-Haq, see Nelsen (2007) p. 116, # 3 ========================

copAMH <-
    (function() { ## to get an environment where  .C  itself is accessible
        C. <- new("acopula", name = "AMH",
                  ## generator
		  psi = function(t,theta) { (1-theta)/(exp(t+0)-theta) },
		  psiInv = function(t,theta) { log((1-theta*(1-t))/t) },
                  ## parameter interval
                  paraInterval = interval("[0,1)"),
                  ## absolute value of generator derivatives
		  psiDabs = function(t, theta, degree=1, n.MC=0, log=FALSE) {
                      if(n.MC > 0) {
                          psiDabsMC(t, family="AMH", theta=theta, degree=degree,
                                    n.MC=n.MC, log=log)
                      } else {
                          if(theta == 0) if(log) return(-t) else return(exp(-t)) # independence
                          ## Note: psiDabs(0, ...) is correct
                          t.et <- theta*exp(-t)
                          if(log)
                              polylog(t.et, s = -degree, method = "neg", log=TRUE)+ log1p(-theta)-log(theta)
                          else
                              polylog(t.et, s = -degree, method = "neg")* (1-theta)/theta
                      }
                  },
                  ## derivatives of the generator inverse
		  psiInvD1abs = function(t, theta, log = FALSE) {
                      if(log) {
                          log1p(-theta)-log(t)-log1p(-theta*(1-t))
                      } else {
                          (1-theta)/(t*(1-theta*(1-t)))
                      }
                  },
                  ## density
		  dacopula = function(u, theta, n.MC=0, log=FALSE) {
                      if(!is.matrix(u)) u <- rbind(u)
                      if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
                      ## f() := NaN outside and on the boundary of the unit hypercube
                      res <- rep.int(NaN, n <- nrow(u))
                      n01 <- apply(u,1,function(x) all(0 < x, x < 1)) # indices for which density has to be evaluated
                      if(!any(n01)) return(res)
                      if(theta == 0) { res[n01] <- if(log) 0 else 1; return(res) } # independence
                      ## auxiliary results
                      u. <- u[n01,, drop=FALSE]
                      u.. <- -theta*(1-u.)
                      sum. <- rowSums(log(u.))
                      sum.. <- rowSums(log1p(u..))
                      ## main part
                      if(n.MC > 0) { # Monte Carlo
                          V <- C.@V0(n.MC, theta)
                          l <- d*log((1-theta)*V) # length = n.MC
                          res[n01] <- colMeans(exp(l + (V-1) %*% t(sum.) - (V+1) %*% t(sum..)))
                          if(log) log(res) else res
                      } else { # explicit
                          Li.arg <- theta*apply(u./(1+u..), 1, prod)
                          Li. <- polylog(Li.arg, s = -d, method = "neg", log=TRUE)
                          res[n01] <- (d+1)*log1p(-theta)-log(theta)+Li.-sum.-sum..
                          if(log) res else exp(res)
                      }
                  },
                  ## score function
                  score = function(u, theta) {
	              if(!is.matrix(u)) u <- rbind(u)
	              if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
                      omu <- 1-u
                      b <- rowSums(omu/(1-theta*omu))
                      h <- theta*apply(u/(1-theta*omu), 1, prod)
                      -(d+1)/(1-theta) - 1/theta + b + (b+1/theta) *
                          polylog(h, s=-(d+1), method="neg") /
                              polylog(h, s=-d, method="neg")
                  },
                  ## nesting constraint
                  nestConstr = function(theta0,theta1) {
                      C.@paraConstr(theta0) &&
                      C.@paraConstr(theta1) && theta1 >= theta0
                  },
                  ## V0 with density dV0 and V01 with density dV01 corresponding to
                  ## LS^{-1}[exp(-V_0psi_0^{-1}(psi_1(t)))]
		  V0 = function(n,theta) rgeom(n, 1-theta) + 1,
		  dV0 = function(x,theta,log = FALSE) dgeom(x-1, 1-theta, log),
                  V01 = function(V0,theta0,theta1) {
                      rnbinom(length(V0),V0,(1-theta1)/(1-theta0))+V0
                  },
                  dV01 = function(x,V0,theta0,theta1,log = FALSE) {
                      stopifnot(length(V0) == 1 || length(x) == length(V0))
                      dnbinom(x-V0,V0,(1-theta1)/(1-theta0),log = log)
                  },
                  ## Kendall's tau
                  tau = tauAMH, ##-> ./aux-acopula.R
                  ## function(th)  1 - 2*((1-th)*(1-th)*log(1-th)+th)/(3*th*th)
                  ## but numerically stable, including, theta -> 0
                  tauInv = function(tau, tol=.Machine$double.eps^0.25, ...) {
                      if(any(tau > 1/3))
                          stop("Impossible for AMH copula to attain a Kendall's tau larger than 1/3")
                      sapply(tau,function(tau) {
                          r <- safeUroot(function(th) tauAMH(th) - tau,
                                         interval = c(0, 1-1e-12),
                                         Sig = +1, tol = tol, check.conv=TRUE, ...)
                          r$root
                      })
                  },
                  ## lower tail dependence coefficient lambda_l
		  lambdaL = function(theta) { 0*theta },
                  lambdaLInv = function(lambda) {
                      if(any(lambda != 0))
                          stop("Any parameter for an Ali-Mikhail-Haq copula gives lambdaL = 0")
                      NA * lambda
                  },
                  ## upper tail dependence coefficient lambda_u
		  lambdaU = function(theta) { 0*theta },
                  lambdaUInv = function(lambda) {
                      if(any(lambda != 0))
                          stop("Any parameter for an Ali-Mikhail-Haq copula gives lambdaU = 0")
                      NA * lambda
                  }
                  )
        C.
    })()# {copAMH}

### ==== Clayton, see Nelsen (2007) p. 116, #1 (slightly simpler form) =========

copClayton <-
    (function() { ## to get an environment where  .C  itself is accessible
        C. <- new("acopula", name = "Clayton",
                  ## generator
		  psi = function(t,theta) { (1+t)^(-1/theta) },
		  psiInv = function(t,theta) { t^(-theta) - 1 },
                  ## parameter interval
                  paraInterval = interval("(0,Inf)"),
                  ## absolute value of generator derivatives
		  psiDabs = function(t, theta, degree=1, n.MC=0, log=FALSE) {
                      if(n.MC > 0) {
                          psiDabsMC(t, family="Clayton", theta=theta, degree=degree,
                                    n.MC=n.MC, log=log)
                      } else {
                          ## Note: psiDabs(0, ...) is correct
                          alpha <- 1/theta
                          res <- lgamma(degree+alpha)-lgamma(alpha)-(degree+alpha)*log1p(t)
                          if(log) res else exp(res)
                      }
                  },
                  ## derivatives of the generator inverse
		  psiInvD1abs = function(t, theta, log = FALSE) {
                      if(log) log(theta)-(1+theta)*log(t) else theta*t^(-(1+theta))
                  },
                  ## density
		  dacopula = function(u, theta, n.MC=0, log=FALSE) {
                      if(!is.matrix(u)) u <- rbind(u)
                      if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
                      ## f() := NaN outside and on the boundary of the unit hypercube
                      res <- rep.int(NaN, n <- nrow(u))
                      n01 <- apply(u,1,function(x) all(0 < x, x < 1)) # indices for which density has to be evaluated
                      if(!any(n01)) return(res)
                      ## auxiliary results
                      u. <- u[n01,, drop=FALSE]
                      l.u <- rowSums(-log(u.))
                      psiI. <- rowSums(C.@psiInv(u., theta))
                      ## main part
                      if(n.MC > 0) { # Monte Carlo
                          l.u.mat <- matrix(rep(l.u, n.MC), nrow=n.MC, byrow=TRUE)
                          V <- C.@V0(n.MC, theta)
                          l <- d*log(theta*V)
                          theta. <- 1 + theta
                          ## stably compute log(colMeans(exp(lx)))
                          lx <- l + theta.*l.u.mat - V %*% t(psiI.) - log(n.MC) # matrix of exponents; dimension n.MC x n ["V x u"]
                          res[n01] <- lsum(lx)
                      } else { # explicit
                          alpha <- 1/theta
                          d.a <- d + alpha
                          res[n01] <- lgamma(d.a)-lgamma(alpha)+ d*log(theta) +
                              (1+theta)*l.u - (d.a)*log1p(psiI.)
                      }
                      if(log) res else exp(res)
                  },
                  ## score function
                  score = function(u, theta) {
	              if(!is.matrix(u)) u <- rbind(u)
	              if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
                      lu <- log(u)
                      alpha <- 1/theta
	              s1 <- (d-sum(1/(theta*(0:(d-1))+1)))*alpha - rowsums(lu)
	              sigma <- rowSums(C.@psiInv(u, theta=theta))
	              s2 <- alpha^2*log(sigma)
	              s3 <- (d+alpha)*rowSums(u^(-theta)*log(u))/sigma
	              s1+s2+s3
                  },
                  ## nesting constraint
                  nestConstr = function(theta0,theta1) {
                      C.@paraConstr(theta0) &&
                      C.@paraConstr(theta1) && theta1 >= theta0
                  },
                  ## V0 with density dV0 and V01 with density dV01 corresponding to
                  ## LS^{-1}[exp(-V_0psi_0^{-1}(psi_1(t)))]
		  V0 = function(n,theta) { rgamma(n, shape = 1/theta) },
		  dV0 = function(x,theta,log = FALSE) dgamma(x, shape = 1/theta, log),
                  V01 = function(V0,theta0,theta1) { retstable(alpha=theta0/theta1, V0) },
                  dV01 = function(x,V0,theta0,theta1,log = FALSE) {
                      stopifnot(length(V0) == 1 || length(x) == length(V0))
                      alpha <- theta0/theta1
                      gamma <- (cos(pi/2*alpha)*V0)^(1/alpha)
                      delta <- V0*(alpha == 1)
		      ## NB: new dstable() is vectorized in (x, gamma, delta) [but not the others]
		      dst <- dstable(x, alpha=alpha, beta = 1, gamma=gamma, delta=delta,
				     pm = 1, log=log, tol = 128* .Machine$double.eps)
                      if(log) V0-x + dst else exp(V0-x) * dst
                  },
                  ## Kendall's tau
		  tau = function(theta) { theta/(theta+2) },
                  tauInv = function(tau) { 2*tau/(1-tau) },
                  ## lower tail dependence coefficient lambda_l
		  lambdaL = function(theta) { 2^(-1/theta) },
                  lambdaLInv = function(lambda) { -1/log2(lambda) },
                  ## upper tail dependence coefficient lambda_u
		  lambdaU = function(theta) { 0*theta },
                  lambdaUInv = function(lambda) {
                      if(any(lambda != 0))
                          stop("Any parameter for a Clayton copula gives lambdaU = 0")
                      NA * lambda
                  }
                  )
        C.
    })()# {copClayton}

### ==== Frank, see Nelsen (2007) p. 116, # 5 ==================================

##' Frank object
copFrank <-
    (function() { ## to get an environment where  .C  itself is accessible
        C. <- new("acopula", name = "Frank",
                  ## generator
		  psi = function(t,theta) {
                      -log1p(expm1(-theta)*exp(0-t))/theta
                      ## == -log(1-(1-exp(-theta))*exp(-t))/theta
                  },
		  psiInv = function(t,theta) {
                      -log(expm1(-theta*t)/expm1(-theta))
                      ## == -log((exp(-theta*t)-1)/(exp(-theta)-1))
                  },
                  ## parameter interval
                  paraInterval = interval("(0,Inf)"),
                  ## absolute value of generator derivatives
		  psiDabs = function(t, theta, degree=1, n.MC=0, log=FALSE) {
                      if(n.MC > 0) {
                          psiDabsMC(t, family="Frank", theta=theta, degree=degree,
                                    n.MC=n.MC, log=log)
                      } else {
                          ## Note: psiDabs(0, ...) is correct
                          p <- -expm1(-theta)
                          if(log)
                              polylog(p*exp(-t), s = -(degree-1), method = "neg", log=TRUE) - log(theta)
                          else
                              polylog(p*exp(-t), s = -(degree-1), method = "neg")/theta
                      }
                  },
                  ## derivatives of the generator inverse
		  psiInvD1abs = function(t, theta, log = FALSE) {
                      if(log) log(theta)-log(expm1(theta*t)) else theta/expm1(theta*t)
                  },
                  ## density
		  dacopula = function(u, theta, n.MC=0, log=FALSE) {
                      if(!is.matrix(u)) u <- rbind(u)
                      if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
                      ## f() := NaN outside and on the boundary of the unit hypercube
                      res <- rep.int(NaN, n <- nrow(u))
                      n01 <- apply(u,1,function(x) all(0 < x, x < 1)) # indices for which density has to be evaluated
                      if(!any(n01)) return(res)
                      ## auxiliary results
                      u. <- u[n01,, drop=FALSE]
                      u.sum <- rowSums(u.)
                      lp <- log1p(-exp(-theta)) # log(1-exp(-theta))
                      lpu <- log1p(-exp(-theta*u.)) # log(1 - exp(-theta * u))
                      lu <- rowSums(lpu)
                      ## main part
                      if(n.MC > 0) { # Monte Carlo
                          V <- C.@V0(n.MC, theta)
                          l <- d*(log(theta*V)-V*lp)
                          rs <- -theta*u.sum
                          lx <- rep(rs, each=n.MC) + l - log(n.MC) + (V-1) %*% t(lu) # (n.MC, nrow(u.))-matrix
                          res[n01] <- lsum(lx)
                          if(log) res else exp(res)
                      } else { # explicit
                          Li.arg <- -expm1(-theta)*exp(rowSums(lpu-lp))
                          Li. <- polylog(Li.arg, s = -(d-1), method = "neg", log=TRUE)
                          res[n01] <- (d-1)*log(theta) + Li. - theta*u.sum - lu
                          if(log) res else exp(res)
                      }
                  },
                  ## score function
                  score = function(u, theta) {
	              if(!is.matrix(u)) u <- rbind(u)
	              if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
                      e <- exp(-theta)
                      e. <- 1-e # 1-e^{-theta}
                      e.u <- exp(-theta*u) # exp(-theta*u)
                      e.u. <- 1 - e.u # 1 - exp(-theta*u)
                      h <- e.*apply(e./e.u, 1, prod)
                      factor <- rowSums(u*e.u/e.u.) - (d-1)*e/e.
                      (d-1)/theta - rowSums(u/e.u.) + factor *
                          polylog(h, s=-d, method="neg") / polylog(h, s=-(d-1), method="neg")
                  },
                  ## nesting constraint
                  nestConstr = function(theta0,theta1) {
                      C.@paraConstr(theta0) &&
                      C.@paraConstr(theta1) && theta1 >= theta0
                  },
                  ## V0 (algorithm of Kemp (1981)) with density dV0 and V01 with density
                  ## dV01 corresponding to LS^{-1}[exp(-V_0psi_0^{-1}(psi_1(t)))]
		  V0 = function(n,theta) rlog(n, -expm1(-theta), exp(-theta)),
		  dV0 = function(x,theta,log = FALSE) {
                      if(any(x != (x <- floor(x + 0.5)))) warning("x must be integer; is rounded with floor(x+0.5) otherwise")
                      if(log) {
                          x*log1p(-exp(-theta))-log(x*theta)
                      } else {
                          p <- 1-exp(-theta)
                          p^x/(x*theta)
                      }
                  },
                  V01 = function(V0, theta0, theta1, rej = 1, approx = 10000) {
                      ## rej method switch: if V0*theta_0*p0^(V0-1) > rej a rejection
                      ## from F01 of Joe is applied (otherwise, the sum is
                      ## sampled via a logarithmic envelope for the summands)
                      ## approx is the largest number of summands before asymptotics is used (see copJoe@V01)
                      ## FIXME: optimal value of rej (for approx = 10000) is not clear yet; rej = 1 is not bad, however
                      ##        lgammacor gives underflow warnings for rej < 1
                      rF01Frank(V0, theta0, theta1, rej, approx)
                  },
                  dV01 = function(x,V0,theta0,theta1,log = FALSE) {
                      stopifnot(length(V0) == 1 || length(x) == length(V0), all(x >= V0))
                      lfactor <- x*log1p(-exp(-theta1))-V0*log1p(-exp(-theta0))
                      res <- lfactor + dsumSibuya(x, V0, theta0/theta1, log=TRUE)
                      if(log) res else exp(res)
                  },
                  ## Kendall's tau; debye_1() is from package 'gsl' :
		  tau = function(theta) {
                      if((l <- length(theta)) == 0) return(numeric(0)) # to work with NULL
                      res <- numeric(l)
                      res[isN <- theta == 0] <- 0 # limiting case
                      res[na <- is.na(theta)] <- NA
                      if(any(i <- !(na | isN)))
                          res[i] <- 1 + 4*(debye_1(theta[i]) - 1)/theta[i]
                      res
                  },
                  tauInv = function(tau, tol = .Machine$double.eps^0.25, ...) {
                      res <- tau
                      res[isN <- res == 0] <- 0 # limiting case
                      res[!isN] <- sapply(res[!isN], function(tau) {
                          r <- safeUroot(function(th) C.@tau(th) - tau,
                                         interval = c(1e-12, 198),
                                         Sig = +1, tol=tol,
                                         check.conv=TRUE, ...)
                          r$root
                      })
                      res
                  },
                  ## lower tail dependence coefficient lambda_l
		  lambdaL = function(theta) { 0*theta },
                  lambdaLInv = function(lambda) {
                      if(any(lambda != 0))
                          stop("Any parameter for a Frank copula gives lambdaL = 0")
                      NA * lambda
                  },
                  ## upper tail dependence coefficient lambda_u
		  lambdaU = function(theta) { 0*theta },
                  lambdaUInv = function(lambda) {
                      if(any(lambda != 0))
                          stop("Any parameter for a Frank copula gives lambdaU = 0")
                      NA * lambda
                  }
                  )
        C.
    })()# {copFrank}

### ==== Gumbel, see Nelsen (2007) p. 116, # 4 =================================

copGumbel <-
    (function() { ## to get an environment where  .C  itself is accessible
        C. <- new("acopula", name = "Gumbel",
                  ## generator
		  psi = function(t,theta) { exp(-t^(1/theta)) },
		  psiInv = function(t,theta) { (-log(t+0))^theta },
                  ## parameter interval
                  paraInterval = interval("[1,Inf)"),
                  ## absolute value of generator derivatives
		  psiDabs = function(t, theta, degree=1, n.MC=0,
                  method=eval(formals(polyG)$method), log = FALSE) {
	              is0 <- t == 0
                      isInf <- is.infinite(t)
	              res <- numeric(n <- length(t))
	              res[is0] <- Inf
                      res[isInf] <- -Inf
                      n0Inf <- !(is0 | isInf)
                      if(all(!n0Inf)) return(if(log) res else exp(res))
                      t. <- t[n0Inf]
                      if(n.MC > 0) {
                          res[n0Inf] <- psiDabsMC(t, family="Gumbel", theta=theta,
                                                  degree=degree, n.MC=n.MC, log=TRUE)
                      } else {
                          if(theta == 1) {
                              res[n0Inf] <- -t. # independence
                          } else {
                              alpha <- 1/theta
                              lt <- log(t.)
                              res[n0Inf] <- -degree*lt -t.^alpha +
                                  polyG(alpha*lt, alpha=alpha, d = degree,
                                        method=method, log=TRUE)
                          }
                      }
                      if(log) res else exp(res)
                  },
                  ## derivatives of the generator inverse
		  psiInvD1abs = function(t, theta, log = FALSE) {
                      if(log) {
                          l.t <- log(t)
                          log(theta)+(theta-1)*log(-l.t)-l.t
                      } else {
                          theta*(-log(t))^(theta-1)/t
                      }
                  },
                  ## density
		  dacopula = function(u, theta, n.MC=0, method=eval(formals(polyG)$method),
                  log = FALSE) {
                      if(!is.matrix(u)) u <- rbind(u)
                      if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
                      ## f() := NaN outside and on the boundary of the unit hypercube
                      res <- rep.int(NaN, n <- nrow(u))
                      n01 <- apply(u,1,function(x) all(0 < x, x < 1)) # indices for which density has to be evaluated
                      if(!any(n01)) return(res)
                      if(theta == 1) { res[n01] <- if(log) 0 else 1; return(res) } # independence
                      ## auxiliary results
                      u. <- u[n01,, drop=FALSE]
                      mlu <- -log(u.) # -log(u)
                      lmlu <- log(mlu) # log(-log(u))
                      psiI. <- rowSums(C.@psiInv(u., theta))
                      ## main part
                      if(n.MC > 0) { # Monte Carlo
                          V <- C.@V0(n.MC, theta)
                          l <- d*log(theta*V)
                          sum. <- rowSums((theta-1)*lmlu + mlu)
                          sum.mat <- matrix(rep(sum., n.MC), nrow=n.MC, byrow=TRUE)
                          ## stably compute log(colMeans(exp(lx)))
                          lx <- l - V %*% t(psiI.) + sum.mat - log(n.MC) # matrix of exponents; dimension n.MC x n ["V x u"]
                          res[n01] <- lsum(lx)
                          if(log) res else exp(res)
                      } else { # explicit
                          alpha <- 1/theta
                          ## compute lx = alpha*log(sum(psiInv(u., theta)))
                          lx <- alpha*log(psiI.)
                          ## ==== former version [start] (numerically slightly more stable but slower) ====
                          ## im <- apply(u., 1, which.max)
                          ## mat.ind <- cbind(seq_len(n), im) # indices that pick out maxima from u.
                          ## mlum <- mlu[mat.ind] # -log(u_max)
                          ## mlum.mat <- matrix(rep(mlum, d), ncol = d)
                          ## lx <- lmlu[mat.ind] + alpha*log(rowSums((mlu/mlum.mat)^theta)) # alpha*log(sum(psiInv(u, theta)))
                          ## ==== former version [end] ====
                          ## compute sum
                          ls. <- polyG(lx, alpha, d, method=method, log=TRUE)-d*lx/alpha
                          ## the rest
                          cop.val <- pnacopula(onacopulaL("Gumbel", list(theta, 1:d)), u.)
                          res[n01] <- d*log(theta) + rowSums((theta-1)*lmlu + mlu) + ls.
                          res[n01] <- if(log) log(cop.val) + res[n01] else cop.val * exp(res[n01])
                          res
                      }
                  },
                  ## score function
                  score = function(u, theta) {
	              if(!is.matrix(u)) u <- rbind(u)
	              if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
                      stop("The score function is currently not implemented for Gumbel copulas")
                  },
                  ## nesting constraint
                  nestConstr = function(theta0,theta1) {
                      C.@paraConstr(theta0) &&
                      C.@paraConstr(theta1) && theta1 >= theta0
                  },
                  ## V0 with density dV0 and V01 with density dV01 corresponding to
                  ## LS^{-1}[exp(-V_0psi_0^{-1}(psi_1(t)))]
		  V0 = function(n,theta) {
                      if(theta == 1) {
                          ## Sample from S(1,1,0,1;1)
                          ## with Laplace-Stieltjes transform exp(-t)
                          rep.int(1., n)
                      } else {
                          alpha <- 1/theta
                          ## Sample from S(alpha,1,(cos(alpha*pi/2))^(1/alpha),0;1)
                          ## with Laplace-Stieltjes transform exp(-t^alpha)
                          rstable1(n, alpha, beta=1,
                                   gamma = cos(alpha*pi/2)^(1/alpha))
                          ## Note: calling sequence:
                          ##       rstable1 -> rstable1C (in rstable1.R) -> rstable1C (in rstable1.R)
                          ##       -> rstable_c (in retstable.c) -> rstable_vec (in retstable.c)
                          ##       -> rstable0 (in retstable.c)
                      }
                  },
		  dV0 = function(x,theta,log = FALSE) C.@dV01(x,1,1,theta,log),
                  V01 = function(V0,theta0,theta1) {
                      alpha <- theta0/theta1
                      if(alpha == 1) {
                          ## Sample from S(1,1,0,V0;1)
                          ## with Laplace-Stieltjes transform exp(-V0*t)
                          V0
                      } else {
                          rstable1(length(V0), alpha, beta=1,
                                   gamma = (cos(alpha*pi/2)*V0)^(1/alpha))
                          ## Sample from S(alpha,1,(cos(alpha*pi/2)V0)^(1/alpha),0;1)
                          ## with Laplace-Stieltjes transform exp(-V0*t^alpha)
                      }
                  },
                  dV01 = function(x,V0,theta0,theta1,log = FALSE) {
                      stopifnot(length(V0) == 1 || length(x) == length(V0))
                      alpha <- theta0/theta1
                      gamma <- (cos(pi/2*alpha)*V0)^(1/alpha)
                      delta <- V0*(alpha == 1)
		      ## NB: new dstable() is vectorized in (x, gamma, delta) [but not the others]
		      dstable(x, alpha=alpha, beta = 1, gamma=gamma, delta=delta, pm = 1, log=log,
			      tol = 128* .Machine$double.eps)
                  },
                  ## Kendall's tau
		  tau = function(theta) { (theta-1)/theta },
                  tauInv = function(tau) { 1/(1-tau) },
                  ## lower tail dependence coefficient lambda_l
		  lambdaL = function(theta) { 0*theta },
                  lambdaLInv = function(lambda) {
                      if(any(lambda != 0))
                          stop("Any parameter for a Gumbel copula gives lambdaL = 0")
                      NA * lambda
                  },
                  ## upper tail dependence coefficient lambda_u
		  lambdaU = function(theta) { 2 - 2^(1/theta) },
                  lambdaUInv = function(lambda) { 1/log2(2-lambda) }
                  )
        C.
    })()# {copGumbel}


### ==== Joe, see Nelsen (2007) p. 116, # 6 ====================================

##' Joe object
copJoe <-
    (function() { ## to get an environment where .C itself is accessible
        C. <- new("acopula", name = "Joe",
                  ## generator
		  psi = function(t,theta) {
                      1 - (-expm1(0-t))^(1/theta)
                      ## == 1 - (1-exp(-t))^(1/theta)
                  },
		  psiInv = function(t,theta) { -log1p(-(1-t)^theta) },
                  ## parameter interval
                  paraInterval = interval("[1,Inf)"),
                  ## absolute value of generator derivatives
		  psiDabs = function(t, theta, degree = 1, n.MC=0,
                  method=eval(formals(polyJ)$method), log = FALSE) {
                      is0 <- t == 0
                      isInf <- is.infinite(t)
                      res <- numeric(n <- length(t))
                      res[is0] <- Inf
                      res[isInf] <- -Inf
                      n0Inf <- !(is0 | isInf)
                      if(all(!n0Inf)) return(if(log) res else exp(res))
                      t. <- t[n0Inf]
                      if(n.MC > 0) {
                          res[n0Inf] <- psiDabsMC(t, family="Joe", theta=theta,
                                                  degree=degree, n.MC=n.MC, log=TRUE)
                      } else {
                          if(theta == 1) {
                              res[n0Inf] <- -t. # independence
                          } else {
                              alpha <- 1/theta
                              mt <- -t.
                              l1mt <- log(-expm1(mt)) # log(1-exp(-t))
                              sum. <- polyJ(mt-l1mt, alpha, degree, method=method, log=TRUE)
                              res[n0Inf] <- -log(theta) + mt - (1-alpha)*l1mt + sum.
                          }
                      }
                      if(log) res else exp(res)
                  },
                  ## derivatives of the generator inverse
		  psiInvD1abs = function(t, theta, log = FALSE) {
                      if(log) {
                          log(theta)+(theta-1)*log1p(-t)-log1p(-(1-t)^theta)
                      } else {
                          theta/((1-t)^(1-theta)-(1-t))
                      }
                  },
                  ## density
		  dacopula = function(u, theta, n.MC=0, method=eval(formals(polyJ)$method),
                  log = FALSE) {
                      if(!is.matrix(u)) u <- rbind(u)
                      if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
                      ## f() := NaN outside and on the boundary of the unit hypercube
                      res <- rep.int(NaN, n <- nrow(u))
                      n01 <- apply(u,1,function(x) all(0 < x, x < 1)) # indices for which density has to be evaluated
                      if(!any(n01)) return(res)
                      if(theta == 1) { res[n01] <- if(log) 0 else 1; return(res) } # independence
                      ## auxiliary results
                      u. <- u[n01,, drop=FALSE]
                      l1_u <- rowSums(log1p(-u.)) # log(1-u)
                      u.. <- (1-u.)^theta # (1-u)^theta
                      lh <- rowSums(log1p(-u..)) # rowSums(log(1-(1-u)^theta)) = log(h)
                      ## main part
                      if(n.MC > 0) { # Monte Carlo
                          V <- C.@V0(n.MC, theta)
                          l <- d*log(theta*V)
                          sum. <- (theta-1)*l1_u
                          sum.mat <- matrix(rep(sum., n.MC), nrow=n.MC, byrow=TRUE)
                          ## stably compute log(colMeans(exp(lx)))
                          lx <- l + (V-1) %*% t(lh) + sum.mat - log(n.MC) # matrix of exponents; dimension n.MC x n ["V x u"]
                          res[n01] <- lsum(lx)
                      } else {
                          alpha <- 1/theta
                          l1_h <- log(-expm1(lh)) # log(1-h)
                          lh_l1_h <- lh - l1_h # log(h/(1-h))
                          res[n01] <- (d-1)*log(theta) + (theta-1)*l1_u -
                              (1-alpha)*l1_h + polyJ(lh_l1_h, alpha, d, method=method,
                                                     log=TRUE)
                      }
                      if(log) res else exp(res)
                  },
                  ## score function
                  score = function(u, theta, method=eval(formals(polyJ)$method)) {
	              if(!is.matrix(u)) u <- rbind(u)
	              if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
                      l1_u <- rowSums(log1p(-u)) # log(1-u)
                      u.th <- (1-u)^theta # (1-u)^theta
                      lh <- rowSums(log1p(-u.th)) # rowSums(log(1-(1-u)^theta)) = log(h)
                      l1_h <- log(-expm1(lh)) # log(1-h)
                      lh_l1_h <- lh - l1_h # log(h/(1-h))
                      b <- rowSums(-l1_u*u.th/(1-u.th))
                      lP <- polyJ(lh_l1_h, alpha, d, method=method, log=TRUE)
                      k <- 1:d
                      alpha <- 1/theta
                      s <- alpha * unlist(lapply(k, function(k.) sum(1/(theta*(1:k.)-1))))
                      ls <- log(s. + (k-1)*exp(log(b) - l1_h))
                      l.a.k <- log(Stirling2.all(d)) + lgamma(k-alpha) - lgamma(1-alpha) + ls
                      lQ <- lsum(l.a.k + (k-1) %*% t(lh_l1_h))
                      (d-1)/theta + rowSums(l1_u) - l1_h/theta^2 + (1-1/theta)*
                          lh_l1_h*b + exp(lQ-lP)
                  },
                  ## nesting constraint
                  nestConstr = function(theta0,theta1) {
                      C.@paraConstr(theta0) &&
                      C.@paraConstr(theta1) && theta1 >= theta0
                  },
                  ## V0 with density dV0 and V01 with density dV01 corresponding to
                  ## LS^{-1}[exp(-V_0psi_0^{-1}(psi_1(t)))]
		  V0 = function(n,theta) rSibuya(n, 1/theta),
		  dV0 = function(x,theta,log = FALSE) {
                      if(log) lchoose(1/theta,x) else abs(choose(1/theta,x))
                  },
                  V01 = function(V0, theta0, theta1, approx = 10000) {
                      ## approx is the largest number of summands before asymptotics is used
                      alpha <- theta0/theta1
                      rF01Joe(V0, alpha, approx)
                  },
		  dV01 =
		  function(x, V0, theta0, theta1,
			   method= eval(formals(dsumSibuya)$method), log = FALSE) {
		      stopifnot(length(V0) == 1 || length(x) == length(V0))
		      ## also holds for theta0 == theta1
		      ## note: this is numerically challenging
		      dsumSibuya(x, V0, theta0/theta1, method=method, log=log)
		  },
                  ## Kendall's tau
                  ## noTerms: even for theta==0, the approximation error is < 10^(-5)
		  tau = function(theta, noTerms=446) {
                      k <- noTerms:1
                      sapply(theta,
                             function(th) {
                                 tk2 <- th*k + 2
                                 1 - 4*sum(1/(k*tk2*(tk2 - th)))
                                 ## ==... (1/(k*(th*k+2)*(th*(k-1)+2)))
                             })
                  },
                  tauInv = function(tau, tol = .Machine$double.eps^0.25, ...) {
                      sapply(tau,function(tau) {
                          r <- safeUroot(function(th) C.@tau(th) - tau,
                                         interval = c(1, 98),
                                         Sig = +1, tol=tol, check.conv=TRUE, ...)
                          r$root
                      })
                  },
                  ## lower tail dependence coefficient lambda_l
		  lambdaL = function(theta) { 0*theta },
                  lambdaLInv = function(lambda) {
                      if(any(lambda != 0))
                          stop("Any parameter for a Joe copula gives lambdaL = 0")
                      NA * lambda
                  },
                  ## upper tail dependence coefficient lambda_u
		  lambdaU = function(theta) { 2-2^(1/theta) },
                  lambdaUInv = function(lambda) { log(2)/log(2-lambda) }
                  )
                  C.
              })()# {copJoe}

### ==== naming stuff ==========================================================

     cNms <- c("copAMH", "copClayton", "copFrank", "copGumbel", "copJoe")
     ## == dput(ls("package:nacopula",pat="^cop"))
     nmsC <- unlist(lapply(cNms, function(.)get(.)@name))
     sNms <- abbreviate(nmsC, 1)
     ## keep these {hidden, for now}:
     c_shortNames <- structure(sNms, names = nmsC)
     c_longNames  <- structure(nmsC, names = sNms)
     c_objNames   <- structure(cNms, names = nmsC)
     rm(cNms, nmsC, sNms)

