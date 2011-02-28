##' <description>
##'
##' @title Outer Power Transformation of an Archimedean Copula
##' @param copbase a "base" copula, i.e. of class "acopula";
##'    must be one of the 5 five predefined families
##' @param thetabase the (univariate) parameter 'theta' for the base copula
##' @return a new "acopula" object; the outer power copula
##' @author Marius Hofert
opower <- function(copbase, thetabase) {
    new("acopula", name = paste("opower", copbase@name, sep=":"),
	## generator
	psi = function(t,theta) { copbase@psi(t^(1/theta), thetabase) },
	psiInv = function(t,theta) { copbase@psiInv(t, thetabase)^theta },
	## parameter interval
	paraInterval = interval("[1,Inf)"),
	## nesting constraint
	nestConstr = function(theta0,theta1) {
	    copbase@paraConstr(theta0) &&
	    copbase@paraConstr(theta1) && theta1 >= theta0
	},
	## V0 and V01
	V0 = function(n,theta) {
	    if(theta == 1) {
		## Sample from S(1,1,0,1;1)
		## with Laplace-Stieltjes transform exp(-t)
		rep.int(1., n)
	    } else {
		V0base <- copbase@V0(n,thetabase) # draw from the base generator
		alpha <- 1/theta
		## Sample from S(alpha,1,(cos(alpha*pi/2))^(1/alpha),0;1)
		## with Laplace-Stieltjes transform exp(-t^alpha)
		S <- rstable1(n, alpha, beta=1,
			      gamma = (cos(alpha*pi/2))^(1/alpha))
		S*V0base^theta
	    }
	},
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
	## Kendall's tau
	tau = function(theta) {
	    1-(1-copbase@tau(thetabase))/theta
	},
	tauInv = function(tau) {
            taubase <- copbase@tau(thetabase)
            if(tau >= taubase) (1-taubase)/(1-tau)
            else {
                stop("The provided tau has to be >= taubase")
                NA * tau
            }
	},
	## lower tail dependence coefficient lambda_l
	lambdaL = function(theta) {
	    if(copbase@name=="Clayton") 2^(-1/(thetabase*theta)) else 0*theta
	},
	lambdaLInv = function(lambda) {
	    if(copbase@name=="Clayton") {
		if(lambda >= 2^(-1/thetabase)) -1/(thetabase*log2(lambda))
		else {
		    stop("The provided lambda has to be >= 2^(-1/thetabase)")
		    NA * lambda
		}
	    } else {
		if(any(lambda != 0))
		    stop("Any parameter for this outer power Archimedean copula gives lambdaL = 0")
		NA * lambda
	    }
	},
	## upper tail dependence coefficient lambda_u
	lambdaU = function(theta) {
	    2 - 2^(1/if(copbase@name %in% c("Gumbel", "Joe")) thetabase*theta else theta)
	},
	lambdaUInv = function(lambda) {
	    if(copbase@name %in% c("Gumbel", "Joe")) {
		if(lambda >= 2-2^(1/thetabase)) 1/(thetabase*log2(2-lambda))
		else {
		    stop("The provided lambda has to be >= 2-2^(1/thetabase)")
		    NA * lambda
		}
	    } else 1/log2(2-lambda)
	})
}
