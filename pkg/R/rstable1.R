##' @param x  numeric vector

##' @return tan(pi *x)	but exact for integer x
tanpi <- function(x) {
    r <- x
    if(any(i <- x == round(x)))
	r[i] <- 0
    io <- which(!i)
    r[io] <- tan(pi* x[io])
    r
}


rstable1 <- function(n, alpha, beta, gamma = 1, delta = 0, pm = 1)
{
    ## Purpose: Random Numbers of Stable distribution in
    ## "0"- or "1" - parametrization
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, based on Diethelm Wuertz's code in fBasics

    ## TODO: Want this to *vectorize* in  alpha, beta, gamma, delta

    stopifnot((la <- length(alpha)) >= 1, (lb <- length(beta)) >= 1,
	      length(gamma) >= 1, length(delta) >= 1,
	      0 < alpha, alpha <= 2, abs(beta) <= 1,
	      length(pm) == 1, pm %in% 0:1)

    p2 <- pi/2
    # Special case (a,b) = (1,0):
    if (all(alpha == 1) && all(beta == 0)) {
	Z <- rcauchy(n)
    }
    else {
	## MM: Nolan(2009) "chapt1.pdf", p.21 has  "Theorem 1.19"
	## -- and attributes that to  'Chambers et al. (1976)'

	## Calculate uniform and exponential distributed random numbers:
	Theta <- pi * (runif(n)-1/2)
	W <- rexp(n)
	##   ^^^^^^ was "-log(runif(n))"; rexp() is faster, giving different numbers

	norm <- alpha != 1 ## TODO:  abs(alpha - 1) > eps.alpha1
	## FIXME(2): ditto for	  | alpha - 1 | << 1

	Z <- numeric(n)
	if(any(norm)) { ## alpha != 1
	    alp <- alpha[norm]; Thet <- Theta[norm]
	    b.tan.pa <- beta[norm]*tanpi(alp/2)
	    th0 <- atan(b.tan.pa) / alp ## == \theta_0
	    ## Now, from Nolan/Chamber's formula, we replace
	    ##	1/(\cos(\alpha\theta_0) \cos\Theta)^{1/\alpha}	with
	    ## c / (\cos\Theta)^{1/\alpha}  where
	    ## c := (1 + \beta*\tan(\pi\alpha/2))^{1/{2\alpha}}
	    ## need to show that c = 1/(\cos(\alpha\theta_0))^{1/\alpha}
	    ## <==> 1 + \beta*\tan(\pi\alpha/2) = 1 / (\cos(\alpha\theta_0))^2
	    c. <- (1 + b.tan.pa^2)^(1/(2*alp))
	    a.tht <- alp*(Thet+th0)
	    Z[norm] <-
		sin(a.tht) * c. / cos(Thet)^(1/alp) *
		    (cos(a.tht-Thet)/W[norm])^((1-alp)/alp) - b.tan.pa
	}
        ## {note that logicals vectorize, indices do *not* so easily}
	if(any(a1 <- !norm)) { ## alpha == 1
	    bet <- beta[a1]; Thet <- Theta[a1]
	    p2 <- pi/2
	    p2.bt <- p2 + bet*Thet
	    Z[a1] <- (p2.bt*tan(Thet) -
		      bet*log((p2*W[a1]*cos(Thet))/p2.bt))/p2
	}
    }

    if(pm == 1) ##  "S" aka "S1" or "1" Parametrization
	delta <- delta + beta * gamma *
	    ifelse(alpha == 1, log(gamma)/p2, tanpi(alpha/2))
    ## else pm == 0 is what we use

    # Result: Location - Scale trafo -- only now using (gamma, delta):
    Z * gamma + delta
}

