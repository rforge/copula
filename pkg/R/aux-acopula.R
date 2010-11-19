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

#### Functions and Methods for "acopula" objects
#### class definition in ./AllClass.R

### ==== Ali-Mikhail-Haq == "AMH" ==============================================

##' @title Ali-Mikhail-Haq ("AMH")'s  tau(theta)
##' @param th
##' @return 1 - 2*((1-th)*(1-th)*log(1-th)+th)/(3*th*th)
##' numerically accurately, notably for th -> 0
##'
##' @author Martin Maechler
tauAMH <- function(th) {
    if(length(th) == 0) return(numeric(0)) # to work with NULL
    r <- th
    lrg <- th > 0.01
    r[lrg] <- (function(t) 1 - 2*((1-t)*(1-t)*log1p(-t) + t)/(3*t*t))(th[lrg])
    if(any(!lrg)) {
	l1 <- !lrg & (ll1 <- th > 2e-4) ## --> k = 7
	r[l1] <- (function(x) 2*x/9*(1+ x*(1/4 + x/10*(1 + x*(1/2 + x/3.5)))))(th[l1])
	l2 <- !ll1 & (ll2 <- th > 1e-5)	 ## --> k = 6
	r[l2] <- (function(x) 2*x/9*(1+ x*(1/4 + x/10*(1 + x/2))))(th[l2])
	l3 <- !ll2 & (ll <- th > 2e-8)	## --> k = 5
	r[l3] <- (function(x) 2*x/9*(1+ x*(1/4 + x/10)))(th[l3])
	irem <- which(!ll)## rest: th <= 2e-8 : k == 4
	r[irem] <- (function(x) 2*x/9*(1+ x/4))(th[irem])
    }
    r
}




### ==== Clayton ===============================================================

##' Note: this is mainly to show that this function can be very well
##' approximated much more simply by just using m <- round(V0).
##' @param V0 numeric vector >= 0
##' @return optimal constant m for the fast rejection algorithm
##' @author Martin Maechler (based on Marius Hofert's code)
m.opt.retst <- function(V0){
    n <- length(V0)
    fV <- floor(V0)
    cV <- ceiling(V0)
    v1 <- fV*exp(V0/fV)
    v2 <- cV*exp(V0/cV)

    m <- integer(n)
    l1 <- (V0 <= 1)
    m[which(l1)] <- 1L

    i2 <- which(!l1) ## those with V0 > 1
    l3 <- (v1[i2] <= v2[i2])
    i3 <- i2[l3]
    m[i3] <- fV[i3]
    i4 <- i2[!l3]
    m[i4] <- cV[i4]
    m
}

##' Sample a random variate St ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)
##' *V_0)^{1/alpha}, V_0*I_{alpha = 1}, I_{alpha != 1}; 1) with
##' Laplace-Stieltjes transform exp(-V_0((1+t)^alpha-1)), see Nolan's book for
##' the parametrization, via an m-fold sum of random variates from
##' \tilde{S}(alpha, 1, (cos(alpha*pi/2)*V_0/m)^{1/alpha}, (V_0/m)
##' *I_{alpha = 1}, I_{alpha != 1}; 1) with Laplace-Stieltjes transform
##' exp(-(V_0/m)*((1+t)^alpha-1)).
##' @param m number of summands, any positive integer
##' @param V0 random variate
##' @param alpha parameter in (0,1]
##' @return St
##' @author Marius Hofert, Martin Maechler
retstablerej <- function(m,V0,alpha){
    gamm. <- (cos(alpha*pi/2)*V0/m)^(1/alpha)
    sum(unlist(lapply(integer(m),
		      function(.) {
			  ## apply standard rejection for sampling
			  ## \tilde{S}(alpha, 1, (cos(alpha*pi/2)
			  ##	*V_0/m)^{1/alpha}, (V_0/m)*I_{alpha = 1},
			  ## h*I_{alpha != 1}; 1) with Laplace-Stieltjes
			  ## transform exp(-(V_0/m)*((h+t)^alpha-h^alpha))
			  repeat {
			      V__ <- rstable1(1, alpha, beta=1, gamma = gamm.)
			      if(runif(1) <= exp(-V__))
				  return(V__)
			  }})
	       ## on acceptance, St_k ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)
	       ## *V_0/m)^{1/alpha}, (V_0/m)*I_{alpha = 1}, h*I_{alpha != 1};
	       ## 1) with Laplace-Stieltjes transform
	       ## exp(-(V_0/m)*((h+t)^alpha-h^alpha))
	       ))
}

##' Sample a vector of random variates St ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)
##' *V_0)^{1/alpha}, V_0*I_{alpha = 1}, h*I_{alpha != 1}; 1) with
##' Laplace-Stieltjes transform exp(-V_0((h+t)^alpha-h^alpha)), see Nolan's book for
##' the parametrization. This procedure calls retstablerej.
##' @param alpha parameter in (0,1]
##' @param V0 vector of random variates
##' @param h non-negative real number
##' @return vector of variates St
##' @author Marius Hofert, Martin Maechler
retstableR <- function(alpha,V0, h = 1){
    n <- length(V0)
    stopifnot(n >= 1, is.numeric(alpha), length(alpha) == 1,
	      0 <= alpha, alpha <= 1) ## <- alpha > 1 ==> cos(pi/2 *alpha) < 0
    ## case alpha==1
    if(alpha == 1) { # alpha == 1 => St corresponds to a point mass at V0 with
	return(V0) # Laplace-Stieltjes transform exp(-V0*t)
    }
    ## else alpha != 1 : call fast rejection algorithm with optimal m
    m <- m.opt.retst(V0)
    mapply(retstablerej, m=m, V0=V0, alpha=alpha)
}

##' Sample random variates St ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)
##' *V_0)^{1/alpha}, V_0*I_{alpha = 1}, h*I_{alpha != 1}; 1) with
##' Laplace-Stieltjes transform exp(-V_0((h+t)^alpha-h^alpha)), see Hofert (2010a).
##' This procedure is more efficient than retstableR since it calls the C
##' function retstable_c.
##' @param alpha parameter in (0,1]
##' @param V0 vector of random variates
##' @param h non-negative real number
##' @param method which method to call ("Marius Hofert", "Luc Devroye")
##' @return vector of variates St
##' @author Martin Maechler
retstableC <- function(alpha, V0, h = 1, method = NULL) {
    n <- length(V0)
    stopifnot(n >= 1, is.numeric(alpha), length(alpha) == 1,
	      0 < alpha, alpha <= 1,
	      is.numeric(h), length(h) == 1, h > 0)
    if(alpha == 1) { # alpha == 1 => St corresponds to a point mass at V0 with
	V0           # Laplace-Stieltjes transform exp(-V0*t)
    }
    else {
	if(is.null(method)) {
	    if(any(diff(V.is.sml <- V0 * h^alpha < 4))) { ## use *both* methods
		r <- numeric(n)
		r[ V.is.sml] <- .Call(retstable_c, V0[ V.is.sml], h = h, alpha, "MH")
		r[!V.is.sml] <- .Call(retstable_c, V0[!V.is.sml], h = h, alpha, "LD")
		return(r)
	    }
	    else
		method <- if(V.is.sml[1]) "MH" else "LD"
	}
	else
	    method <- match.arg(method, c("MH","LD"))
	.Call(retstable_c, V0, h = h, alpha, method)
    }
}

## switch to make fast C code the default
if(FALSE)
    retstable <- retstableR
retstable <- retstableC

### ==== Frank =================================================================

##' Random number generator for a Log(p) distribution, R version
##' @param n number of random variates to be generated
##' @param p parameter in (0,1)
##' @return vector of random variates from Log(p)
##' @author Marius Hofert, Martin Maechler
rlogR <- function(n,p) {
    stopifnot((n <- as.integer(n)) >= 0, 0 < p, p < 1)
    vec <- numeric(n)
    if(n >= 1) {
	u <- runif(n)
	l1 <- u > p
	vec[l1] <- 1
	i2 <- which( !l1 ) # of shorter length, say n2
	q2 <- 1-(1-p)^runif(length(i2)) # length n2
	l3 <- u[i2] < q2*q2
	i3 <- i2[l3]
	vec[i3] <- floor(1+abs(log(u[i3])/log(q2[l3])))
	l4 <- u[i2] > q2
	vec[i2[l4]] <- 1
	l5 <- ! (l3 | l4) # q2^2 <= u[i2] <= q2
	vec[i2[l5]] <- 2
    }
    vec
}

##' Random number generator for a Log(p) distribution, C version
##' @param n number of random variates to be generated
##' @param p parameter in (0,1)
##' @return vector of random variates from Log(p)
##' @author Martin Maechler
rlog <- function(n,p) {
    stopifnot(n >= 0,  0 < p, p < 1)
    .Call(rLog_c, n, p)
}

##' Sample V ~ F with Laplace-Stieltjes transform
##' (1-(1-exp(-t)*(1-e^(-theta1)))^alpha)/(1-e^(-theta0))
##' via the algorithm of Hofert (2010a). R version.
##' @param p parameter 1-e^(-theta1)
##' @param alpha parameter theta0/theta1 in (0,1]
##' @param theta0_le_1 in {0,1} with 1 if and only if theta0 <= 1
##' @return V
##' @author Marius Hofert, Martin Maechler
rejFFrankR <- function(p,alpha,theta0_le_1) {
    if(theta0_le_1) {
	repeat{
	    U <- runif(1)
	    X <- rlog(1,p)
	    if(U*(X-alpha) <= 1/beta(X,1-alpha)) break
	}
    } else {
	repeat{
	    U <- runif(1)
	    X <- rFJoe(1,alpha)
	    if(U <= p^(X-1)) break
	}
    }
    X
}

##' Vectorize rejFFrankR. Generate a vector of variates V ~ F with
##' Laplace-Stieltjes transform (1-(1-exp(-t)*(1-e^(-theta1)))^alpha)
##' /(1-e^(-theta0)). R version.
##' @param n length of the vector of random variates
##' @param theta0 parameter theta0 in (0,infinity)
##' @param theta1 parameter theta1 in [theta0, infinity)
##' @return vector of random variates from F
##' @author Martin Maechler
rFFrankR <- function(n,theta0,theta1) {
    sapply(rep.int(-expm1(-theta1), n), rejFFrankR,
	   alpha = theta0/theta1,
	   theta0_le_1 = (theta0 <= 1))
}

##' FIXME: provide title... + maybe put in C (if it works)
##' @param V0 one realization of V0
##' @param theta0 parameter theta0 in (0,infinity)
##' @param theta1 parameter theta1 in [theta0, infinity)
##' @return vector of random variates from F01
##' @author Marius Hofert
rF01FrankR <- function(V0,theta0,theta1) {
	p1 <- -expm1(-theta1)
    	repeat{
	    U <- runif(1)
	    X <- copJoe@V01(V0,theta0,theta1) # sample V01 of Joe
	    if(U <= p1^X) break
	}
	X
}

##' Generate a vector of variates V ~ F with Laplace-Stieltjes transform
##' (1-(1-exp(-t)*(1-e^(-theta1)))^alpha)/(1-e^(-theta0)). C version.
##' @param n length of the vector of random variates
##' @param theta0 parameter theta0 in (0,infinity)
##' @param theta1 parameter theta1 in [theta0, infinity)
##' @return vector of random variates from F
##' @author Martin Maechler
rFFrank <- function(n,theta0,theta1) {
    .Call(rFFrank_c, n, theta0, theta1);
}

### ==== Joe ===================================================================

##' Sample V ~ F with F(n) = 1-1/(n*B(n,1-alpha)), n in IN, with
##' Laplace-Stieltjes transform 1-(1-exp(-t))^alpha via the algorithm of
##' Hofert (2010a). R version.
##' @param n  sample size
##' @param alpha parameter
##' @return vector of random variates V
##' @author Marius Hofert, Martin Maechler
rFJoeR <- function(n,alpha) {
    stopifnot((n <- as.integer(n)) >= 0)
    V <- numeric(n)
    if(n >= 1) {
	if(alpha == 1) {
	    V[] <- 1
	} else {
	    u <- runif(n)
	    ## FIXME(MM): (for alpha not too close to 1): re-express using 1-u
	    l1 <- u <= alpha
	    V[l1] <- 1
	    i2 <- which(!l1)
	    Ginv <- ((1-u[i2])*gamma(1-alpha))^(-1/alpha)
	    floorGinv <- floor(Ginv)
	    l3 <- (1-1/(floorGinv*beta(floorGinv,1-alpha)) < u[i2])
	    V[i2[l3]] <- ceiling(Ginv[l3])
	    i4 <- which(!l3)
	    V[i2[i4]] <- floorGinv[i4]
	}
    }
    V
}

##' Sample V ~ F with F(n) = 1-1/(n*B(n,1-alpha)), n in IN, with
##' Laplace-Stieltjes transform 1-(1-exp(-t))^alpha via the algorithm of
##' Hofert (2010a). C version.
##' @param n  sample size
##' @param alpha parameter
##' @return vector of random variates V
##' @author Martin Maechler
rFJoe <- function(n,alpha) {
    stopifnot(is.numeric(n), n >= 0)
    .Call(rFJoe_c, n, alpha)
}

### ==== other numeric utilities ===============================================

##' Compute Stirling numbers of the 1st kind
##'
##' s(n,k) = (-1)^{n-k} times
##' the number of permutations of 1,2,…,n with exactly k cycles
##'
##' NIST DLMF 26.8 --> http://dlmf.nist.gov/26.8
##' @title  Stirling Numbers of the 1st Kind
##' @param n
##' @param k
##' @return s(n,k)
##' @author Martin Maechler
Stirling1 <- function(n,k)
{
    ## NOTA BENE: There's no "direct" method available here
    stopifnot(length(n) == 1, length(k) == 1)
    if (k < 0 || n < k) stop("'k' must be in 0..n !")
    if(n == 0) return(1)
    if(k == 0) return(0)
    S1 <- function(n,k) {
	if(k == 0 || n < k) return(0)
	if(is.na(S <- St[[n]][k])) {
	    ## s(n,k) = s(n-1,k-1) - (n-1) * s(n-1,k) for all n, k >= 0
	    St[[n]][k] <<- S <- if(n1 <- n-1L)
                S1(n1, k-1) - n1* S1(n1, k) else 1
	}
	S
    }
    if(compute <- (nt <- length(St <- get("S1.tab", .nacopEnv))) < n) {
	## extend the "table":
	length(St) <- n
	for(i in (nt+1L):n) St[[i]] <- rep.int(NA_real_, i)
    }
    else compute <- is.na(S <- St[[n]][k])
    if(compute) {
	S <- S1(n,k)
	## store it back:
	assign("S1.tab", St, envir = .nacopEnv)
    }
    S
}

##' Full Vector of Stirling Numbers of the 1st Kind
##' @title  Stirling1(n,k) for all k = 1..n
##' @param n
##' @return the same as sapply(1:n, Stirling1, n=n)
##' @author Martin Maechler
Stirling1.all <- function(n)
{
    stopifnot(length(n) == 1)
    if(!n) return(numeric(0))
    if(get("S1.full.n", .nacopEnv) < n) {
	assign("S1.full.n", n, envir = .nacopEnv)
	unlist(lapply(seq_len(n), Stirling1, n=n))
    }
    else get("S1.tab", .nacopEnv)[[n]]
}

##' Compute Stirling numbers of the 2nd kind
##'
##' S^{(k)}_n = number of ways of partitioning a set of $n$ elements into $k$
##'	non-empty subsets
##' (Abramowitz/Stegun: 24,1,4 (p. 824-5 ; Table 24.4, p.835)
##'   Closed Form : p.824 "C."
##' @title  Stirling Numbers of the 2nd Kind
##' @param n
##' @param k
##' @param method
##' @return S(n,k) = S^{(k)}_n
##' @author Martin Maechler, "direct": May 28 1992
Stirling2 <- function(n,k, method = c("lookup.or.store","direct"))
{
    stopifnot(length(n) == 1, length(k) == 1)
    if (k < 0 || n < k) stop("'k' must be in 0..n !")
    switch(match.arg(method),
	   "direct" = {
	       sig <- rep(c(1,-1)*(-1)^k, length= k+1) # 1 for k=0; -1 1 (k=1)
	       k <- 0:k # (!)
	       ga <- gamma(k+1)
	       round(sum( sig * k^n /(ga * rev(ga))))
	   },
	   "lookup.or.store" = {
               if(n == 0) return(1) ## else:
               if(k == 0) return(0)
	       S2 <- function(n,k) {
		   if(k == 0 || n < k) return(0)
		   if(is.na(S <- St[[n]][k]))
		       ## S(n,k) = S(n-1,k-1) + k * S(n-1,k)   for all n, k >= 0
		       St[[n]][k] <<- S <- if(n1 <- n-1L)
                           S2(n1, k-1) + k* S2(n1, k) else 1 ## n = k = 1
		   S
	       }
	       if(compute <- (nt <- length(St <- get("S2.tab", .nacopEnv))) < n) {
		   ## extend the "table":
		   length(St) <- n
		   for(i in (nt+1L):n) St[[i]] <- rep.int(NA_real_, i)
	       }
	       else compute <- is.na(S <- St[[n]][k])
	       if(compute) {
		   S <- S2(n,k)
		   ## store it back:
		   assign("S2.tab", St, envir = .nacopEnv)
	       }
	       S
	   })
}

##' Full Vector of Stirling Numbers of the 2nd Kind
##' @title  Stirling2(n,k) for all k = 1..n
##' @param n
##' @return the same as sapply(1:n, Stirling2, n=n)
##' @author Martin Maechler
Stirling2.all <- function(n)
{
    stopifnot(length(n) == 1)
    if(!n) return(numeric(0))
    if(get("S2.full.n", .nacopEnv) < n) {
	assign("S2.full.n", n, envir = .nacopEnv)
	unlist(lapply(seq_len(n), Stirling2, n=n))
    }
    else get("S2.tab", .nacopEnv)[[n]]
}

## Our environment for tables etc:  no hash, as it will contain *few* objects:
.nacopEnv <- new.env(parent=emptyenv(), hash=FALSE)
assign("S2.tab", list(), envir = .nacopEnv) ## S2.tab[[n]][k] == S(n, k)
assign("S1.tab", list(), envir = .nacopEnv) ## S1.tab[[n]][k] == s(n, k)
assign("S2.full.n", 0  , envir = .nacopEnv)
assign("S1.full.n", 0  , envir = .nacopEnv)


##' From   http://en.wikipedia.org/wiki/Polylogarithm
##' 1. For integer values of the polylogarithm order, the following
##'   explicit expressions are obtained by repeated application of z·∂/∂z
##'   to Li1(z):
##'
##'     {Li}_{1}(z) = -\ln(1-z)
##'     {Li}_{0}(z) = {z \over 1-z}
##'     {Li}_{-1}(z) = {z \over (1-z)^2}
##'     {Li}_{-2}(z) = {z \,(1+z) \over (1-z)^3}
##'     {Li}_{-3}(z) = {z \,(1+4z+z^2) \over (1-z)^4}
##'     {Li}_{-4}(z) = {z \,(1+z) (1+10z+z^2) \over (1-z)^5}.
##'
##' Accordingly the polylogarithm reduces to a ratio of polynomials in
##' z, and is therefore a rational function of z, for all nonpositive
##' integer orders. The general case may be expressed as a finite sum:
##'
##' {Li}_{-n}(z) = \left(z \,{\partial \over \partial z} \right)^n \frac{z}{1-z}=
##'     = \sum_{k=0}^n k! \,S(n+1,k+1) \left({z \over {1-z}} \right)^{k+1}
##' \ \ (n=0,1,2,\ldots),
##'
##' where S(n,k) are the Stirling numbers of the second
##' kind. Equivalent formulae applicable to negative integer orders are
##' (Wood 1992, § 6):
##'
##'  {Li}_{-n}(z) = (-1)^{n+1} \sum_{k=0}^n k! \,S(n+1,k+1) \left({{-1} \over {1-z}} \right)^{k+1} \
##'     (n=1,2,3,\ldots),
##'
##' Compute the polylogarithm function \eqn{Li_s(z)}
##'
##' @title Polylogarithm Li_s(z)
##' @param z numeric or complex vector
##' @param s complex number; current implementation is aimed at s \in 0,-1,..
##' @param method a string specifying the algorithm to be used
##' @param logarithm
##' @return numeric/complex vector as \code{z}
##' @author Martin Maechler
polylog <- function(z,s, method = c("sum","negint-s_Stirling"), logarithm=FALSE,
		    ## for "sum" -- this is more for experiments etc:
		    n.sum)
{
    if((nz <- length(z)) == 0 || (ns <- length(s)) == 0)
	return((z+s)[FALSE])# of length 0
    if(logarithm)
	stop("'logarithm = TRUE' is not yet implemented")
    stopifnot(length(s) == 1) # for now
    method <- match.arg(method)
    switch(method,
	   "sum" = {
	       stopifnot((Mz <- Mod(z)) <= 1, Mz < 1 | Re(s) > 1,
			 n.sum > 99, length(n.sum) == 1)
	       n <- n.sum:1L # reversed in order, so we sum *small* terms first
	       if(nz == 1)
		   sum(z^n / n^s)
	       else
		   colSums(outer(n, z, function(n,z) z^n) / rep.int(n, nz)^s)
	   },
	   "negint-s_Stirling" = {
	       stopifnot(s == as.integer(s), s <= 1)
	       if(s == 1) return(-log1p(-z)) ## -ln(1 -z)
	       r <- z/(1 - z)
	       ## if(s == 0) return(r)
	       n <- -as.integer(s)
	       k1 <- seq_len(n+1)# == k+1, k = 0...n
	       fac.k <- cumprod(c(1, seq_len(n)))
               S.n1.k1 <- Stirling2.all(n+1) ## == Stirling2(n+1, k+1)
	       if(nz == 1)
		   sum(	   fac.k * S.n1.k1  *  r^k1 )
	       else
		   colSums(fac.k * S.n1.k1 * outer(k1, r, function(k,r) r^k))
	   }, stop("invalid 'method':", method))
}

### ==== other NON-numerics ====================================================

##' Function which computes psiDAbs via Monte Carlo
##' @param t evaluation points
##' @param family Archimedean family
##' @param theta parameter value
##' @param degree order of derivative
##' @param MC Monte Carlo sample size
##' @param log if TRUE the log of psiDAbs is returned
psiDAbsMC <- function(t, family, theta, degree = 1, MC, log = FALSE){
	V0. <- getAcop(family)@V0(MC,theta)
        l.V0. <- degree*log(V0.)
        summands <- function(t) mean(exp(-V0.*t + l.V0.))
        res <- unlist(lapply(t,summands))
        if(log) log(res) else res
}

##' Function for setting the parameter in an acopula
##' @param x acopula
##' @param value parameter value
##' @param na.ok logical indicating if NA values are ok for theta
##' @return acopula with theta set to value
##' @author Martin Maechler
setTheta <- function(x, value, na.ok = TRUE) {
    stopifnot(is(x, "acopula"),
	      is.numeric(value) | (ina <- is.na(value)))
    if(ina) {
	if(!na.ok) stop("NA value, but 'na.ok' is not TRUE")
	value <- NA_real_
    }
    if(ina || x@paraConstr(value)) ## parameter constraints are fulfilled
	x@theta <- value
    else
	stop("theta (=", format(value), ")  does not fulfil paraConstr()")
    x
}


##' Construct "paraConstr" function from an "interval"
##' @param int interval
##' @return parameter constraint function
##' @author Martin Maechler
mkParaConstr <- function(int){
    stopifnot(is(int, "interval")) # for now
    is.o <- int@open
    eL <- substitute(LL <= theta, list(LL = int[1])); if(is.o[1]) eL[[1]] <-
	as.symbol("<")
    eR <- substitute(theta <= RR, list(RR = int[2])); if(is.o[2]) eR[[1]] <-
	as.symbol("<")
    bod <- substitute(length(theta) == 1 && LEFT && RIGHT,
		      list(LEFT = eL, RIGHT= eR))
    as.function(c(alist(theta=), bod), parent.env(environment()))
    ## which is a fast version of
    ## r <- function(theta) {}
    ## environment(r) <- parent.env(environment())
    ## body(r) <- bod
    ## r
}

printAcopula <- function(x, slots = TRUE, indent = 0,
			 digits = getOption("digits"), width = getOption("width"), ...)
{
    cl <- class(x)
    cld <- getClassDef(cl)
    stopifnot(indent >= 0, extends(cld, "acopula"))
    ch.thet <- {
	if(!all(is.na(x@theta)))## show theta
	    paste(", theta= (",
		  paste(sapply(x@theta, format, digits=digits), collapse=", "),
		  ")", sep="")
	else ""
    }
    bl <- paste(rep.int(" ",indent), collapse="")
    cat(sprintf('%sArchimedean copula ("%s"), family "%s"%s\n',
		bl, cl, x@name, ch.thet))
    if(slots) {
	nms <- slotNames(cld)
	nms <- nms[!(nms %in% c("name", "theta"))]
	i2 <- indent+2
	cat(bl, " It contains further slots, named\n",
	    paste(strwrap(paste(dQuote(nms),collapse=", "),
			  width = 0.95 * (width-2), indent=i2, exdent=i2),
		  collapse="\n"), "\n",
	    sep="")
    }
    invisible(x)
}
setMethod(show, "acopula", function(object) printAcopula(object))

## This is now exported & has help file --> ../man/printNacopula.Rd :
printNacopula <-
    function(x, labelKids = NA, deltaInd = if(identical(labelKids,FALSE)) 5 else 3,
             indent.str="",
	     digits = getOption("digits"), width = getOption("width"), ...)
{
    cl <- class(x)
    stopifnot(deltaInd >= 0, is.character(indent.str), length(indent.str) == 1,
              extends(cl, "nacopula"))
    mkBlanks <- function(n) paste(rep.int(" ", n), collapse="")
    bl <- mkBlanks(nIS <- nchar(indent.str))

    ## cat(sprintf(" __deltaInd = %d, nIS = %d__ ", deltaInd, nIS))
    ch1 <- sprintf("%sNested Archimedean copula (\"%s\"), with ",
                   indent.str, cl)
    ch2 <- if(length(c.j <- x@comp)) {
	sprintf("slot \n%s'comp'   = %s", bl,
		paste("(",paste(c.j, collapse=", "),")", sep=""))
    } else "empty slot 'comp'"
    cat(ch1, ch2, sprintf("  and root\n%s'copula' = ", bl), sep="")
    printAcopula(x@copula, slots=FALSE, digits=digits, width=width, ...)
    nk <- length(kids <- x@childCops)
    if(nk) {
	cat(sprintf("%sand %d child copula%s\n", bl, nk, if(nk > 1)"s" else ""))
	doLab <- if(is.na(labelKids)) nk > 1 else as.logical(labelKids)
	paste0 <- function(...) paste(..., sep="")
	if(doLab) {
	    hasNms <- !is.null(nms <- names(kids))
	    lab <- if(hasNms) paste0(nms,": ") else paste0(seq_len(nk),") ")
	}
        bl <- mkBlanks(nIS + deltaInd)
	for(ic in seq_along(kids))
	    printNacopula(kids[[ic]], deltaInd=deltaInd,
			  indent.str = paste0(bl, if(doLab) lab[ic]),
			  labelKids=labelKids, digits=digits, width=width, ...)
    }
    else
	cat(sprintf("%sand *no* child copulas\n", bl))
    invisible(x)
}

setMethod(show, "nacopula", function(object) printNacopula(object))


##' @title Get one of our "acopula" family objects by name
##' @param family either character string (short or longer form of
##'	 copula family name) or an "acopula" family object
##' @param check logical indicating if the class of the return value should
##' be checked.
##' @return one of our "acopula" objects
##' @author Martin Maechler
getAcop <- function(family, check=TRUE) {
    if(is.character(family)) {
        stopifnot(length(family) == 1)
        if(nchar(family) <= 2)          # it's a short name
            family <- c_longNames[family]
        COP <- get(c_objNames[family])  # envir = "package:nacopula"
        if(check && !is(COP, "acopula"))
            stop(paste("invalid acopula-family object, family=",family))
        COP
    } else if(is(family, "acopula"))
        family
    else stop("'family' must be an \"acopula\" object or family name")
}
