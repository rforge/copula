sinc <- function(x) .Call(sinc_c, x)

A..Z <- function(x, alpha, I.alpha = 1 - alpha)
    .Call(nacopula:::A__c, x, alpha, I.alpha)

##' <description>
##'
##' @title tan(pi *x) exact for integer x
##' @param x  numeric vector
##' @return   numeric vector as x
##' @author Martin Maechler
tanpi <- function(x) {
    r <- x
    if(any(i <- x == round(x)))
	r[i] <- 0
    io <- which(!i)
    r[io] <- tan(pi* x[io])
    r
}

##' @title  cos(pi/2 * x)
##' @param x  numeric vector

##' @return cos(pi/2 *x)	but exact for integer x
##' @author Martin Maechler
cospi2 <- function(x) {
    r <- x
    if(any(i <- x == round(x))) {
        xi.4 <- x[i] %% 4
	r[xi.4 == 0     ] <-  1
	r[xi.4 == 2     ] <- -1
	r[xi.4 %% 2 != 0] <-  0
    }
    io <- which(!i)
    r[io] <- cos(pi/2 * x[io])
    r
}


rstable1R <- function(n, alpha, beta, gamma = 1, delta = 0, pm = 1)
{
    ## Purpose: Random Numbers of Stable distribution in
    ## "0"- or "1" - parametrization
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, based on Diethelm Wuertz's code in fBasics

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

        a.is.vec <- (la > 1) # alpha is "vector" (not scalar)
        if(a.is.vec) {
            ## if 'alpha' is not scalar, make sure that lengths of
            ##		alpha, beta, Theta, W  are all equal (== n)
            alpha <- rep(alpha, length.out = n)
            beta  <- rep(beta,  length.out = n)
        }

	norm <- alpha != 1 ## TODO:  abs(alpha - 1) > eps.alpha1
	## FIXME(2): ditto for	  | alpha - 1 | << 1

	Z <- numeric(n)
	if(any(norm)) { ## alpha != 1
	    alp <- alpha[norm]; Thet <- Theta[norm]
	    b.tan.pa <- beta[norm]*tanpi(alp/2)
	    th0 <- atan(b.tan.pa) / alp ## == \theta_0
	    ## Now, from Nolan/Chambers' formula, we replace
	    ##	1/(\cos(\alpha\theta_0) \cos\Theta)^{1/\alpha}	with
	    ## c / (\cos\Theta)^{1/\alpha}  where
	    ## c := (1 + (\beta*\tan(\pi\alpha/2))^2)^{1/{2\alpha}}
	    ## need to show that c = 1/(\cos(\alpha\theta_0))^{1/\alpha}
	    ## <==> 1 + (\beta*\tan(\pi\alpha/2))^2 = 1 / (\cos(\alpha\theta_0))^2
            ##  and that's true,  as  1 + (tan(y))^2 = 1 / (cos(y))^2  for
            ##   y = \alpha\theta_0 = \arc\tan(\beta*\tan(\pi\alpha/2))
	    c. <- (1 + b.tan.pa^2)^(1/(2*alp))
	    a.tht <- alp*(Thet+th0)
	    Z[norm] <-
		sin(a.tht) * c. / cos(Thet)^(1/alp) *
		    (cos(a.tht-Thet)/W[norm])^((1-alp)/alp)
	}
        ## {note that logicals vectorize, indices do *not* so easily}
	if(any(a1 <- !norm)) { ## alpha == 1
	    bet <- beta[a1]; Thet <- Theta[a1]
	    p2 <- pi/2
	    p2.bt <- p2 + bet*Thet
	    Z[a1] <- (p2.bt*tan(Thet) - bet*log((p2*W[a1]*cos(Thet))/p2.bt))/p2
	}
    }

    if(pm == 0)
        ## delta_1 := delta_0 - .. [Nolan, chapt.1, (1.7)], since above 1-parametr.
	delta <- delta - gamma * {
	    if(a.is.vec) {
		d.of <- numeric(n)
		d.of[norm] <- b.tan.pa
		d.of[a1	 ] <- beta * log(gamma)/p2
		d.of
	    }
	    else { ## alpha is scalar
		if(norm) b.tan.pa else beta * log(gamma)/p2
	    }
	}

    # Result: Location - Scale trafo -- only now using (gamma, delta):
    Z * gamma + delta
}

rstable1C <- function(n, alpha, beta, gamma = 1, delta = 0, pm = 1)
{
    ## Purpose: Random Numbers of Stable distribution in
    ## Using the  'pm == 1'  parametrization only
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler

    stopifnot((la <- length(alpha)) >= 1, (lb <- length(beta)) >= 1,
	      length(gamma) >= 1, length(delta) >= 1,
	      0 < alpha, alpha <= 2, abs(beta) <= 1,
              length(pm) == 1, pm %in% 0:1)

    if(beta == 1 && pm == 1) ## use the C version
        .Call(rstable_c, n, alpha) * gamma + delta
    else rstable1R(n, alpha=alpha, beta=beta,
                   gamma=gamma, delta=delta, pm=pm)
}

rstable1 <- rstable1R

### ?? BUG ???
## --- I can test compare the above and all looks fine,
##    see ks.test() in ../tests/rstable-ex.R
##
## if(FALSE) ## but if we do
rstable1 <- rstable1C
## the correlation (tau) tests fail for "Gumbel"
##>> ../tests/NAC-experi.R
##>> stopifnot(max(abs(corCheckGumbel[["cor"]]-trCorr)) < eps.tau)

