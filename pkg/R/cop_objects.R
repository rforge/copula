###' List of supported Archimedean copulas

##' FIXME: Not "nice" that the names of the copula objects must be be used
##'	inside some of their components.
##' Possible solution:
##'	use *same* environment for  (tau, tauInv, paraConstr, nestConstr)
##' NOTA BENE:  Write psi(), tau(), ... functions such that they *vectorize*
##' ---------   *and* do work even for (non-numeric) NULL argument
##'          	{now checked in "ACopula" validityMethod -- see ./AllClass.R }

##' ==== Ali-Mikhail-Haq, see Nelsen (2007) p. 116, # 3 ========================

copAMH <-
    new("ACopula", name = "AMH",
        ##' generator
        psi = function(t,theta) { (1-theta)/(exp(t+0)-theta) },
        psiInv = function(t,theta) { log((1-theta*(1-t))/t) },
        ##' parameter interval
        paraInterval = interval("[0,1)"),
        ##' nesting constraint
        nestConstr = function(theta0,theta1) {
            copAMH@paraConstr(theta0) &&
            copAMH@paraConstr(theta1) && theta1 >= theta0
        },
        ##' V0 and V01
        V0 = function(n,theta) { rgeom(n, 1-theta) + 1 },
        V01 = function(V0,theta0,theta1) {
            rnbinom(length(V0),V0,(1-theta1)/(1-theta0))+V0
        },
        ##' Kendall's tau
        tau = function(theta) {
            1 - 2*((1-theta)*(1-theta)*log(1-theta)+theta)/(3*theta*theta)
        },
        tauInv = function(tau, tol = .Machine$double.eps^0.25, ...) {
            if(any(tau > 1/3))
                stop("Impossible for AMH copula to attain a Kendall's tau larger
 than 1/3")
            sapply(tau,function(tau) {
		r <- safeUroot(function(th) copAMH@tau(th) - tau,
			       interval = c(1e-12, 1-1e-12), Sig = +1, 
                               tol = tol, ...)
                r$root ##' FIXME: check for convergence
            })
        },
        ##' lower tail dependence coefficient lambda_l
        lTDC = function(theta) { 0*theta },
        lTDCInv = function(lambda) {
            if(any(lambda != 0))
                stop("Any parameter for an Ali-Mikhail-Haq copula gives zero 
lower tail dependence coefficient")
            NA * lambda
        },
        ##' upper tail dependence coefficient lambda_u
        uTDC = function(theta) { 0*theta },
        uTDCInv = function(lambda) {
            if(any(lambda != 0))
                stop("Any parameter for an Ali-Mikhail-Haq copula gives zero 
upper tail dependence coefficient")
            NA * lambda
        })

##' ==== Clayton, see Nelsen (2007) p. 116, #1 (slightly simpler form here) ====

copClayton <-
    new("ACopula", name = "Clayton",
        ##' generator
        psi = function(t,theta) { (1+t)^(-1/theta) },
        psiInv = function(t,theta) { t^(-theta) - 1 },
        ##' parameter interval
        paraInterval = interval("(0,Inf)"),
        ##' nesting constraint
        nestConstr = function(theta0,theta1) {
            copClayton@paraConstr(theta0) &&
            copClayton@paraConstr(theta1) && theta1 >= theta0
        },
        ##' V0 and V01
        V0 = function(n,theta) { rgamma(n, shape = 1/theta) },
        V01 = function(V0,theta0,theta1) { retstable(alpha=theta0/theta1, V0) },
        ##' Kendall's tau
        tau = function(theta) { theta/(theta+2) },
        tauInv = function(tau) { 2*tau/(1-tau) },
        ##' lower tail dependence coefficient lambda_l
        lTDC = function(theta) { 2^(-1/theta) },
        lTDCInv = function(lambda) { -1/log2(lambda) },
        ##' upper tail dependence coefficient lambda_u
        uTDC = function(theta) { 0*theta },
        uTDCInv = function(lambda) {
            if(any(lambda != 0))
                stop("Any parameter for a CLayton copula gives zero upper tail 
dependence coefficient")
            NA * lambda
        })

##' ==== Frank, see Nelsen (2007) p. 116, # 5 ==================================

##' Frank object
copFrank <-
    new("ACopula", name = "Frank",
        ##' generator
	psi = function(t,theta) {
	    -log1p(expm1(-theta)*exp(0-t))/theta
	    ##' == -log(1-(1-exp(-theta))*exp(-t))/theta
	},
	psiInv = function(t,theta) {
	    -log(expm1(-theta*t)/expm1(-theta))
	    ##' == -log((exp(-theta*t)-1)/(exp(-theta)-1))
	},
        ##' parameter interval
        paraInterval = interval("(0,Inf)"),
        ##' nesting constraint
        nestConstr = function(theta0,theta1) {
            copFrank@paraConstr(theta0) &&
            copFrank@paraConstr(theta1) && theta1 >= theta0
        },
        ##' V0 (algorithm of Kemp (1981)) and V01
        V0 = function(n,theta) { rlog(n,1-exp(-theta)) },
	V01 = function(V0,theta0,theta1) {
            ##' FIXME: how to approximate when V0 large? not yet solved 
            ##' theoretically
            sapply(lapply(V0, rFFrank, theta0=theta0,theta1=theta1), sum)
	},
        ##' Kendall's tau; debye_1() is from package 'gsl' :
        tau = function(theta) 1 + 4*(debye_1(theta) - 1)/theta,
	tauInv = function(tau, tol = .Machine$double.eps^0.25, ...) {
	    sapply(tau, function(tau) {
		r <- safeUroot(function(th) copFrank@tau(th) - tau,
			       interval = c(0.001,100), Sig = +1, tol=tol, ...)
                ##' FIXME: check for convergence
                r$root
            })
        },
        ##' lower tail dependence coefficient lambda_l
        lTDC = function(theta) { 0*theta },
        lTDCInv = function(lambda) {
            if(any(lambda != 0))
                stop("Any parameter for a Frank copula gives zero lower tail 
dependence coefficient")
            NA * lambda
        },
        ##' upper tail dependence coefficient lambda_u
        uTDC = function(theta) { 0*theta },
        uTDCInv = function(lambda) {
            if(any(lambda != 0))
                stop("Any parameter for a Frank copula gives zero upper tail 
dependence coefficient")
            NA * lambda
        })

##' ==== Gumbel, see Nelsen (2007) p. 116, # 4 =================================

copGumbel <-
    new("ACopula", name = "Gumbel",
        ##' generator 
        psi = function(t,theta) { exp(-t^(1/theta)) },
        psiInv = function(t,theta) { (-log(t+0))^theta },
        ##' parameter interval
        paraInterval = interval("[1,Inf)"),
        ##' nesting constraint
        nestConstr = function(theta0,theta1) {
            copGumbel@paraConstr(theta0) &&
            copGumbel@paraConstr(theta1) && theta1 >= theta0
        },
        ##' V0 and V01
	V0 = function(n,theta) {
	    if(theta == 1) {
		##' Sample from S(1,1,0,1;1)
		##' with Laplace-Stieltjes transform exp(-t)
		rep.int(1., n)
	    } else {
		alpha <- 1/theta
		##' Sample from S(alpha,1,(cos(alpha*pi/2))^(1/alpha),0;1)
		##' with Laplace-Stieltjes transform exp(-t^alpha)
		rstable1(n, alpha, beta=1,
			 gamma = cos(alpha*pi/2)^(1/alpha))
	    }
	},
	V01 = function(V0,theta0,theta1) {
	    alpha <- theta0/theta1
	    if(alpha == 1) {
		##' Sample from S(1,1,0,V0; 1)
		##' with Laplace-Stieltjes transform exp(-V0*t)
		V0
	    } else {
		rstable1(length(V0), alpha, beta=1,
			 gamma = (cos(alpha*pi/2)*V0)^(1/alpha))
		##' Sample from S(alpha,1,(cos(alpha*pi/2)V0)^(1/alpha),0;1)
		##' with Laplace-Stieltjes transform exp(-V0*t^alpha)
	    }
	},
        ##' Kendall's tau
        tau = function(theta) { (theta-1)/theta },
        tauInv = function(tau) { 1/(1-tau) },
        ##' lower tail dependence coefficient lambda_l
        lTDC = function(theta) { 0*theta },
        lTDCInv = function(lambda) {
            if(any(lambda != 0))
                stop("Any parameter for a Gumbel copula gives zero lower tail 
dependence coefficient")
            NA * lambda
        },
        ##' upper tail dependence coefficient lambda_u
        uTDC = function(theta) { 2 - 2^(1/theta) },
        uTDCInv = function(lambda) { 1/log2(2-lambda) }
        ) #' {copGumbel}

##' ==== Joe, see Nelsen (2007) p. 116, # 6 ====================================

##' Joe object
copJoe <-
    new("ACopula", name = "Joe",
        ##' generator
        psi = function(t,theta) {
            1 - (-expm1(0-t))^(1/theta)
            ##' == 1 - (1-exp(-t))^(1/theta)
        },
        psiInv = function(t,theta) { -log1p(-(1-t)^theta) },
        ##' parameter interval
        paraInterval = interval("[1,Inf)"),
        ##' nesting constraint
        nestConstr = function(theta0,theta1) {
            copJoe@paraConstr(theta0) &&
            copJoe@paraConstr(theta1) && theta1 >= theta0
        },
        ##' V0 and V01
        V0 = function(n,theta) rFJoe(n, 1/theta),
        V01 = function(V0,theta0,theta1, approx=100000) {
            alpha <- theta0/theta1
            ia <- (V0 > approx)
            ie <- !ia
            V01 <- V0 #' numeric(length(V0))
            V01[ia] <- V0[ia]^(1/alpha) *
                rstable1(sum(ia),alpha, beta=1, 
                         gamma = cos(alpha*pi/2)^(1/alpha),
                         delta = as.integer(alpha==1))
            V01[ie] <- sapply(lapply(V0[ie], rFJoe, alpha=alpha),sum)
            V01
        },
        ##' Kendall's tau
        ##' noTerms: even for theta==0, the approximation error is < 10^(-5)
        tau = function(theta, noTerms=446) {
            k <- seq_len(noTerms)
            sapply(theta,
                   function(th) {
                       tk2 <- th*k + 2
                       1 - 4*sum(1/(k*tk2*(tk2 - th)))
                       ##' ==... (1/(k*(th*k+2)*(th*(k-1)+2)))
                   })
        },
        tauInv = function(tau, tol = .Machine$double.eps^0.25, ...) {
            sapply(tau,function(tau) {
		r <- safeUroot(function(th) copJoe@tau(th) - tau,
			       interval = c(1.001, 100), Sig = +1, tol=tol, ...)
                ##' FIXME: check for convergence
                r$root
            })
        },
        ##' lower tail dependence coefficient lambda_l
        lTDC = function(theta) { 0*theta },
        lTDCInv = function(lambda) {
            if(any(lambda != 0))
                stop("Any parameter for a Joe copula gives zero lower tail
 dependence coefficient")
            NA * lambda
        },
        ##' upper tail dependence coefficient lambda_u
        uTDC = function(theta) { 2-2^(1/theta) },
        uTDCInv = function(lambda) { log(2)/log(2-lambda) }
        ) ##' {copJoe}

##' ==== other stuff ===========================================================

cNms <- c("copAMH", "copClayton", "copFrank", "copGumbel", "copJoe")
##' == dput(ls("package:nacopula",pat="^cop"))
nmsC <- unlist(lapply(cNms, function(.)get(.)@name))
sNms <- abbreviate(nmsC, 1)
##' keep these {hidden, for now}:
c_shortNames <- structure(sNms, names = nmsC)
c_longNames  <- structure(nmsC, names = sNms)
c_objNames   <- structure(cNms, names = nmsC)
rm(cNms, nmsC, sNms)
