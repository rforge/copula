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
##' @author Martin Maechler
tauAMH <- function(th) {
    if(length(th) == 0) return(numeric(0)) # to work with NULL
    r <- th
    na <- is.na(th)
    lrg <- (th > 0.01) & !na
    r[lrg] <- (function(t) 1 - 2*((1-t)*(1-t)*log1p(-t) + t)/(3*t*t))(th[lrg])
    if(any(!lrg & !na)) {
	l1 <- !lrg & !na & (ll1 <- th > 2e-4) ## --> k = 7
	r[l1] <- (function(x) 2*x/9*(1+ x*(1/4 + x/10*(1 + x*(1/2 + x/3.5)))))(th[l1])
	l2 <- !ll1 & !na & (ll2 <- th > 1e-5)	 ## --> k = 6
	r[l2] <- (function(x) 2*x/9*(1+ x*(1/4 + x/10*(1 + x/2))))(th[l2])
	l3 <- !ll2 & !na & (ll <- th > 2e-8)	## --> k = 5
	r[l3] <- (function(x) 2*x/9*(1+ x*(1/4 + x/10)))(th[l3])
	irem <- which(!ll & !na)## rest: th <= 2e-8 : k == 4
	r[irem] <- (function(x) 2*x/9*(1+ x/4))(th[irem])
    }
    r
}


### ==== Clayton ===============================================================

##' Note: this is mainly to show that this function can be very well
##' approximated much more simply by just using m <- round(V0).
##'
##' @title Optimal constant for fast rejection
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

### ==== fast rejection for fixed m, R version ====

##' Sample a random variate St ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)
##' *V_0)^{1/alpha}, V_0*I_{alpha = 1}, I_{alpha != 1}; 1) with
##' Laplace-Stieltjes transform exp(-V_0((1+t)^alpha-1)), see Nolan's book for
##' the parametrization, via an m-fold sum of random variates from
##' \tilde{S}(alpha, 1, (cos(alpha*pi/2)*V_0/m)^{1/alpha}, (V_0/m)
##' *I_{alpha = 1}, I_{alpha != 1}; 1) with Laplace-Stieltjes transform
##' exp(-(V_0/m)*((1+t)^alpha-1)). This is a building block for the fast rejection.
##'
##' @title Sample an exponentially tilted stable distribution as an m-fold sum
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

### ==== fast rejection, R version ====

##' Sample a vector of random variates St ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)
##' *V_0)^{1/alpha}, V_0*I_{alpha = 1}, h*I_{alpha != 1}; 1) with
##' Laplace-Stieltjes transform exp(-V_0((h+t)^alpha-h^alpha)), see Nolan's book for
##' the parametrization, with the fast rejection. This procedure calls retstablerej.
##'
##' @title Sampling an exponentially tilted stable distribution
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

### ==== state-of-the-art: fast rejection + Luc's algorithm, C version ====

##' Sample random variates St ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)
##' *V_0)^{1/alpha}, V_0*I_{alpha = 1}, h*I_{alpha != 1}; 1) with
##' Laplace-Stieltjes transform exp(-V_0((h+t)^alpha-h^alpha)), see Nolan's book for
##' the parametrization, with the fast rejection.
##' This procedure is more efficient than retstableR since it calls the C
##' function retstable_c and uses both the fast rejection and Luc Devroye's algorithm.
##'
##' @title Efficiently sampling an exponentially tilted stable distribution
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
retstable <- retstableC # retstable is by default retstableC


### ==== Frank =================================================================

### ==== sampling a logarithmic distribution, R version ====

##' Random number generator for a Log(p) distribution with the algorithm "LK" of
##' Kemp (1981), R version.
##'
##' @title Sample a Log(p) distribution
##' @param n number of random variates to be generated
##' @param p parameter in (0,1)
##' @param Ip = 1 - p_ (possibly more accurate) -- use, if p ~= 1
##' @return vector of random variates from Log(p)
##' @author Marius Hofert, Martin Maechler
rlogR <- function(n, p, Ip = 1-p) {
    if(missing(p)) p <- 1 - Ip
    stopifnot((n <- as.integer(n)) >= 0,
              0 < p, p <= 1, 0 < Ip)
    vec <- numeric(n)
    if(n >= 1) {
	u <- runif(n)
	l1 <- u > p
	vec[l1] <- 1
	i2 <- which( !l1 ) # of shorter length, say n2
	q2 <- 1-(Iq2 <- Ip^runif(length(i2))) # length n2
	l3 <- u[i2] < q2*q2
	i3 <- i2[l3]
	vec[i3] <- floor(1+abs(log(u[i3])/log1p(-Iq2[l3])))
	l4 <- u[i2] > q2
	vec[i2[l4]] <- 1
	l5 <- ! (l3 | l4) # q2^2 <= u[i2] <= q2
	vec[i2[l5]] <- 2
    }
    vec
}

### ==== state-of-the art: sampling a logarithmic distribution, C version ====

##' Random number generator for a Log(p) distribution with the algorithm "LK" of
##' Kemp (1981), C version.
##'
##' @title Efficiently sampling a Log(p) distribution
##' @param n number of random variates to be generated
##' @param p parameter in (0,1)
##' @param Ip = 1 - p_ (possibly more accurate)
##' @return vector of random variates from Log(p)
##' @author Martin Maechler
rlog <- function(n, p, Ip = 1-p) {
    if(missing(p)) p <- 1 - Ip
    stopifnot(n >= 0, 0 < p, p <= 1, 0 < Ip)
    .Call(rLog_vec_c, n, p, Ip)
}

### ==== state-of-the art: sampling F01Frank, C version ====

##' Generate a vector of variates V01 ~ F01 with Laplace-Stieltjes transform
##' ((1-(1-exp(-t)*(1-e^(-theta1)))^alpha)/(1-e^(-theta0)))^V0.
##'
##' @title Efficiently sampling V01 for Frank
##' @param V0 vector of random variates from F0
##' @param theta0 parameter theta0 in (0,infinity)
##' @param theta1 parameter theta1 in [theta0, infinity)
##' @param rej method switch: if V0*theta_0*p0^(V0-1) > rej a rejection
##'        from F01 of Joe is applied (otherwise, the sum is
##'        sampled via a logarithmic envelope for the summands)
##' @param approx largest number of summands before asymptotics is used
##' @return vector of random variates V01
##' @author Marius Hofert
rF01Frank <- function(V0, theta0, theta1, rej, approx) {
    .Call(rF01Frank_vec_c, V0, theta0, theta1, rej, approx)
}

### ==== wrapper for inner distribution F for Frank ====

##' Generate a vector of variates V ~ F with Laplace-Stieltjes transform
##' (1-(1-exp(-t)*(1-e^(-theta1)))^alpha)/(1-e^(-theta0)).
##'
##' @title Sampling F for Frank
##' @param n number of variates from F
##' @param theta0 parameter theta0 in (0,infinity)
##' @param theta1 parameter theta1 in [theta0, infinity)
##' @param rej method switch: if theta_0 > rej a rejection
##'        from Joe (Sibuya distribution) is applied (otherwise, a logarithmic
##' 	   envelope is used)
##' @return vector of random variates V
##' @author Marius Hofert
rFFrank <- function(n, theta0, theta1, rej)
    rF01Frank(rep(1,n),theta0,theta1,rej,1) # approx = 1 (dummy)


### ==== Gumbel ================================================================

### ==== compute the coefficients for polyG ====

##' Compute the coefficients a_{dk}(theta) involved in the generator derivatives
##' and the copula density of a Gumbel copula
##'
##' @title Coefficients of the polynomial involved in the generator derivatives
##'        and density for Gumbel
##' @param d number of coefficients, d >= 1
##' @param alpha parameter (1/theta)
##' @param method method slot:
##'        binomial.coeff: currently best method available, uses binomial coefficients,
##'                        only critical for large dependencies
##'        sort:           compute the coefficients via exp(log()), pulling out the maximum, and sort
##'        horner:         uses polynEval
##'        direct:         brute force approach
##'        dV01.Joe:       uses dV01.Joe
##' @param log boolean which determines if the logarithm is returned
##' @param verbose logical for method == sort
##' @return a_{dk}(theta) = (-1)^{d-k}\sum_{j=k}^d alpha^j * s(d,j) * S(j,k)
##' FIXME: serious numerical problems, e.g., for (d = 100, alpha = 0.8)
##' @author Marius Hofert und Martin Maechler
coeffG <- function(d, alpha,
		   method = c("sort", "horner", "direct", "dV01.Joe"),
		   log = FALSE, verbose = FALSE)
{
    stopifnot(d >= 1)
    a <- numeric(d) # for the a_{dk}(theta)'s
    method <- match.arg(method)
    switch(method,
	   "sort" = {
	       ls <- log(abs(Stirling1.all(d))) # log(|s(d, i)|), i=1:d
	       lS <- lapply(1:d, function(n) log(Stirling2.all(n)))
	       ##-> lS[[n]][i] contains log(S(n,i)), i = 1,...,n
               wrong.sign <- integer()
	       for(k in 1:d) { # deal with a_{dk}(theta)
		   j <- k:d
		   ## compute b_j's, j = k,..,d
		   b <- j * log(alpha) + ls[j] +
		       unlist(lapply(j, function(i) lS[[i]][k]))
		   b.max <- max(b) # max{b_j}
		   exponents <- b - b.max # exponents
		   ## compute critical sum (known to be positive)
		   exps <- exp(exponents) # (d-k+1)-many
		   even <- if(k == d) NULL else seq(2, d-k+1, by=2)
		   odd <- seq(1, d-k+1, by=2)
		   sum.neg <- sum(sort(exps[even]))
		   sum.pos <- sum(sort(exps[odd]))
		   sum. <- sum.pos - sum.neg
		   a[k] <- if(log) b.max + log(sum.) else exp(b.max)*sum.
		   if(sum.neg > sum.pos) {
		       if(verbose) message("sum.neg > sum.pos for k = ", k)
		       wrong.sign <- c(wrong.sign, k)
		   }
	       }
	       if(length(wrong.sign))
		   attr(a, "wrong.signs") <- wrong.sign
	       a
	   },
	   "dV01.Joe" = {
	       ## coefficients via dV01 of Joe
	       ## a_{dk}(theta) = d!/k! * copJoe@dV01(d, k, theta0=1, theta1=1/alpha)
	       k <- 1:d
	       ck <- ## c_k := d!/k!
		   if(log) c(0,cumsum(log(d:2)))[d:1]
		   else c(1,cumprod(d:2))[d:1]
	       ## copJoe @ dV01() is already vectorized (somewhat):
	       p <- copJoe@dV01(rep(d,d), k, theta0 = 1, theta1 = 1/alpha,
				log=log) ## FIXME: call dJoe
	       if(log) p + ck else p * ck
	   },
	   "horner" = {
	       s.abs <- abs(Stirling1.all(d))
	       S <- lapply(1:d, Stirling2.all)
               ## S[[n]][i] contains S(n,i), i = 1,...,n
	       k <- 1:d
	       pol <- vapply(k, function(k.) {
		   j <- 0:(d-k.)
		   ## compute coefficients (c_k = |s(d,j+k)|*S(j+k,k))
		   c.j <- s.abs[j+k.] *
		       unlist(lapply(j, function(i) S[[i+k.]][k.]))
		   polynEval(c.j, -alpha)
	       }, NA_real_)

	       if(log) k*log(alpha) + log(pol) else alpha^k * pol
	   },
	   "direct" = {
	       s <- Stirling1.all(d) # s(d,1), ..., s(d,d)
	       k <- 1:d
	       S <- lapply(k, Stirling2.all) # S[[k]][n] contains S(k,n), n = 1,...,k
	       vapply(k, function(k.) {
		   j <- k.:d
		   ## extract a column of Stirling2 numbers:
		   S. <- unlist(lapply(j, function(i) S[[i]][k.]))
		   sm <- sum(alpha^j * s[j]*S.)
		   if(log) log(abs(sm)) else (-1)^(d-k.)*sm
	       }, NA_real_)
	   },
       {stop(sprintf("unsupported method '%s' in coeffG", method))}
           ) ## switch()
} ## coeffG()

### ==== compute the polynomial for Gumbel ====

##' Compute the polynomial involved in the generator derivatives and the
##' copula density of a Gumbel copula
##'
##' @title Polynomial involved in the generator derivatives and density for Gumbel
##' @param lx = log(x); where x: evaluation point (vector);
##'        e.g., for copGumbel@dacopula, lx = alpha*log(rowSums(psiInv(u)))
##'        where u = (u_1,..,u_d) is the evaluation point of the density of Joe's copula)
##' @param alpha parameter (1/theta)
##' @param d number of summands, >= 1
##' @param method a string, one of
##'   pois:           uses ppois
##'   binomial.coeff: uses binomial coefficients, only critical for large dependencies
##'   sort:           compute the coefficients via exp(log()),
##' 			  pulling out the maximum, and sort
##'   horner:         uses polynEval
##'   direct:         brute force approach
##'   dV01.Joe:       uses dV01.Joe
##' @param log boolean which determines if the logarithm is returned
##' @return \sum_{k=1}^d  a_{dk}(\theta)  x ^ k
##'       = \sum_{k=1}^d  a_{dk} *     exp(lx*k)
##'  where a_{dk}(theta)
##'       = (-1)^{d-k}\sum_{j=k}^d \theta^{-j} s(d,j) S(j,k)
##'       = (d!/k!)\sum_{l=1}^k (-1)^{d-l} \binom{k}{l}\binom{\alpha l}{d}
##' @author Marius Hofert
polyG <- function(lx, alpha, d, method=c("pois", "binomial.coeff", "sort",
				"horner", "direct", "dV01.Joe"), log=FALSE)
{
    k <- 1:d
    method <- match.arg(method)
    switch(method,
           "pois" =
       { ### FIXME: improve speed
           ## determine signs of the falling factorials
           signs <- (-1)^k* (2*(floor(alpha*k) %% 2) - 1)

           ## build list of b's
           n <- length(lx)
           x <- exp(lx) ## e^lx = x
           lppois <- outer(d-k, x, FUN=ppois, log.p=TRUE) # a (d x n)-matrix; log(ppois(d-k, x))
           llx <- outer(k, lx) # also a (d x n)-matrix; k*lx
           labsPoch <- unlist(lapply(k, function(j) sum(log(abs(alpha*j-0:(d-1))) ) )) # log|(alpha*k)_d|
           lfac <- lfactorial(k)
           ## build matrix of exponents
           B <- llx + lppois + rep(labsPoch - lfac, n) + rep(x, each = d)
           max.B <- apply(B, 2, max)
           ## pull out maximum and sum the rest
           res <- max.B + log(as.vector(signs %*% exp(B - rep(max.B, each=d))))
           if(log) res else exp(res)
       },
           "binomial.coeff" =
       { ### FIXME: improve speed --- careful!! - see end of ../demo/logL-vis.R

           ## determine signs of the falling factorials
           ## s <- unlist(lapply(alpha*k, function(z) prod(z-(0:(d-1)))))
           ## signs <- (-1)^(d-k) * sign(s) ## see  ../misc/sign-polyG.R
           ## signs  <- (-1)^(d-k) * (-1)^d * (2*(floor(alpha*k) %% 2) - 1)
           ## signs  =  (-1)^k              * (2*(floor(alpha*k) %% 2) - 1)
           signs <- (-1)^k* (2*(floor(alpha*k) %% 2) - 1)

           ## build list of b's
           n <- length(lx)
           x <- exp(lx) ## e^lx = x
           k1 <- k-1L # = 0:(d-1)
           lppois <- outer(d-k, x, FUN=ppois, log.p=TRUE) # a (d x n)-matrix; log(ppois(d-k, x))
           llx <- outer(k, lx)		# also a (d x n)-matrix; k*lx
           labsPoch <- vapply(k, function(j) sum(log(abs(alpha*j-k1))),
                              NA_real_)# log|(alpha*k)_d|
           lfac <- lfactorial(k)
           ## build matrix of exponents
           B <- llx + lppois + rep(labsPoch - lfac, n) + rep(x, each = d)
           max.B <- apply(B, 2, max)
           ## pull out maximum and sum the rest
           res <- max.B + log(as.vector(signs %*% exp(B - rep(max.B, each=d))))
           if(log) res else exp(res)
       },

           "sort" =, "horner" =, "direct" =, "dV01.Joe" =
       {
           ## note: these methods are all know to show numerical deficiencies
           if(d > 220) stop("d > 220 not yet supported") # would need Stirling2.all(d, log=TRUE)
           ## compute the log of the coefficients:
           a.dk <- coeffG(d, alpha, method=method)
           l.a.dk <- log(a.dk) # note: theoretically, a.dk > 0 but due to numerical issues, this might not always be the case
           ## evaluate the sum
           ## for this, create a matrix B with (k,i)-th entry
           ## B[k,i] = log(a_{dk}(theta)) + k * lx[i],
           ##          where k in {1,..,d}, i in {1,..,n} [n = length(lx)]
           B <- l.a.dk + outer(k, lx)
           if(log){
               ## compute log(colSums(exp(B))) stably (no overflow) with the idea of
               ## pulling out the maxima
               max.B <- apply(B, 2, max)
               max.B + log(colSums(exp(B - rep(max.B, each = d))))
           }else colSums(exp(B))
       },
           stop(sprintf("unsupported method '%s' in polyG",
                        method))
           ) # end{switch}
}


### ==== Joe ===================================================================

### ==== sampling a Sibuya(alpha) distribution, R version ====

##' Sample V from a Sibuya(alpha) distribution with cdf F(n) = 1-1/(n*B(n,1-alpha)),
##' n in IN, with Laplace-Stieltjes transform 1-(1-exp(-t))^alpha via the
##' algorithm of Hofert (2011), Proposition 3.2. R version.
##'
##' @title Sampling Sibuya(alpha) distributions
##' @param n  sample size
##' @param alpha parameter
##' @return vector of random variates V
##' @author Marius Hofert, Martin Maechler
rSibuyaR <- function(n,alpha) {
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

### ==== state-of-the-art: sampling a Sibuya(alpha) distribution, C version ====

##' Sample V from a Sibuya(alpha) distribution with cdf F(n) = 1-1/(n*B(n,1-alpha)),
##' n in IN, with Laplace-Stieltjes transform 1-(1-exp(-t))^alpha via the
##' algorithm of Hofert (2011), Proposition 3.2. C version.
##'
##' @title Efficiently sampling Sibuya(alpha) distributions
##' @param n  sample size
##' @param alpha parameter
##' @return vector of random variates V
##' @author Martin Maechler
rSibuya <- function(n,alpha) {
    stopifnot(is.numeric(n), n >= 0)
    .Call(rSibuya_vec_c, n, alpha)
}

### ==== state-of-the art: sampling F01Joe, C version ====

##' Generate a vector of variates V01 ~ F01 with Laplace-Stieltjes transform
##' ((1-(1-exp(-t))^alpha))^V0. Bridge to R. Used, e.g., to draw several variates
##' from rF01Joe.
##'
##' @title Sampling F01 for Joe's family
##' @param V0 vector of random variates from F0
##' @param parameter alpha = theta0/theta1 in (0,1]
##' @param approx largest number of summands before asymptotics is used
##' @return vector of random variates V01
##' @author Marius Hofert
rF01Joe <- function(V0, alpha, approx) {
    .Call(rF01Joe_vec_c, V0, alpha, approx)
}

### ==== wrapper for inner distribution F for Joe ====

##' Generate a vector of variates V ~ F with Laplace-Stieltjes transform
##' 1-(1-exp(-t))^alpha.
##'
##' @title Sampling F for Joe
##' @param n number of variates from F
##' @param parameter alpha = theta0/theta1 in (0,1]
##' @return vector of random variates V
##' @author Marius Hofert
rFJoe <- function(n, alpha) rSibuya(n, alpha)

### ==== polynomial evaluation for Joe ====

##' Inner probability mass function for a nested Joe copula
##'
##' @title Inner probability mass function for a nested Joe copula
##' @param x vector (or number) of natural numbers
##' @param k vector (or number) of natural numbers
##' @param alpha parameter in (0,1]
##' @param method method applied
##'        log:      proper log computation based on lssum
##'        direct:   brute-force evaluation of the sum and its log
##'        exp.log:  similar to method = "log", but without *proper/intelligent* log
##' @param log boolean which determines if the logarithm is returned
##' @return \sum_{j=1}^k choose(k,j)*choose(alpha*j,x)*(-1)^(x-j)
##' @author Marius Hofert
##' note: - this is a probability mass function in x, where x in {k, k+1, ...}
##'       - numerically challenging, e.g., dJoe(100, 96, 0.01) < 0 for all methods
dJoe <- function(x, k, alpha,
		 method=c("log", "direct", "exp.log"), log=FALSE)
{
    if(alpha == 1) outer(x, k, FUN="==") # if alpha = 1 then the result is x == k
    l.x <- length(x)
    l.k <- length(k)
    res <- matrix(, nrow=l.x, ncol=l.k)
    ## k1
    k1 <- k==1
    ## numerically more stable computation for k=1 :
    res[,k1] <- if(log) lchoose(alpha, x) else abs(choose(alpha, x))
    ## !k1
    if(any(!k1)) {
	method <- match.arg(method)
	switch(method,
               "log" = {
	           ## computes *proper* log based on lssum [stops if numerical problem appears]
                   ## determine the matrix of signs of (alpha*j,x)*(-1)^(x-j), j in {1,..,m}
                   ## note: this does not depend on x!
                   m <- max(k[!k1])
                   signs <- unlist(lapply(1:m, function(j){
                       z <- alpha*j
                       if(z == floor(z)) 0 else (-1)^(j-ceiling(z))
                   }))
                   ## for one pair of x and k:
                   one.arg <- function(z){ # z = (x,k[!k1])
	               if(z[1] < z[2]) return(-Inf) # = log(0)
                       j <- 1:z[2]
                       lxabs <- lchoose(z[2], j) + lchoose(alpha*j, z[1]) # build incredient for lssum
                       lssum(lxabs, signs[1:z[2]]) # call lssum
                   }
                   sum. <- matrix(apply(expand.grid(x, k[!k1]), 1, FUN=one.arg), nrow=l.x)
                   res[,!k1] <- if(log) sum. else exp(sum.)
               },
               "direct" = {
	           ## brute force evaluation of the sum and its log
                   one.arg <- function(z){ # z = (x,k[!k1])
	               if(z[1] < z[2]) return(0)
                       j <- 1:z[2]
                       sum(choose(z[2],j)*choose(alpha*j,z[1])*(-1)^(z[1]-j))
                   }
                   sum. <- matrix(apply(expand.grid(x, k[!k1]), 1, FUN=one.arg), nrow=l.x)
                   res[,!k1] <- if(log) log(sum.) else sum.
               },
               "exp.log" = {
	           ## similar to method = "log", but without *proper/intelligent* log
                   ## and inefficient due to the signs (old version)
                   one.arg <- function(z){ # z = (x,k[!k1])
                       if(z[1] < z[2]) return(0)
                       j <- 1:z[2] # indices of the summands
                       signs <- (-1)^(j+z[1])
                       ## determine the signs of choose(j*alpha,z[1]) for each component of j
                       to.subtract <- 0:(z[1]-1)
                       signs.choose <- unlist(lapply(j,function(l){
                           prod(sign(l*alpha-to.subtract))}
                                                     ))
                       signs <- signs*signs.choose
                       binom.coeffs <- exp(lchoose(z[2],j)+lchoose(j*alpha,z[1]))
                       sum(signs*binom.coeffs)
                   }
                   sum. <- matrix(apply(expand.grid(x, k[!k1]), 1, FUN=one.arg), nrow=l.x)
                   res[,!k1] <- if(log) log(sum.) else sum.
               },
           {stop(sprintf("unsupported method '%s' in dJoe", method))})
    }
    if(l.x == 1 || l.k == 1) as.vector(res) else res
}

### ==== polynomial evaluation for Joe ====

##' Compute the polynomial involved in the generator derivatives and the
##' copula density of a Joe copula
##'
##' @title Polynomial involved in the generator derivatives and density for Joe
##' @param lx (log) evaluation point (lx is meant to be log(x) for some x which
##'        was used earlier; e.g., for copJoe@dacopula, lx = log(h(u)/(1-h(u))) for
##'        h(u) = \prod_{j=1}^d(1-(1-u_j)^theta), where u = (u_1,..,u_d) is the
##'        evaluation point of the density of Joe's copula)
##' @param alpha parameter (1/theta)
##' @param d number of summands
##' @param method different methods, can be
##'        "log.poly" intelligent log version
##'        "log1p"    additonally uses log1p
##'        "poly"     brute force log version
##' @param log boolean which determines if the logarithm is returned
##' @return \sum_{k=1}^d a_{dk}(theta) exp((k-1)*lx),
##'         where a_{dk}(theta) = S(d,k)*(k-1-alpha)_{k-1} = S(d,k)*Gamma((1:d)-alpha)/Gamma(1-alpha)
##' @author Marius Hofert and Martin Maechler
polyJ <- function(lx, alpha, d, method=c("log.poly","log1p","poly"), log=FALSE){
    ## compute the log of the coefficients a_{dk}(theta)
    if(d > 220) stop("d > 220 not yet supported")# would need Stirling2.all(d, log=TRUE)
    k <- 1:d
    l.a.k <- log(Stirling2.all(d)) + lgamma(k-alpha) - lgamma(1-alpha) # log(a_{dk}(theta)), k = 1,..,d
    ## FIXME: maybe (!) use Horner (see polyG)
    ## evaluate polynomial via exp( log(<poly>) )
    ## for this, create a matrix B with (k,i)-th entry B[k,i] = log(a_{dk}(theta)) + (k-1) * lx[i],
    ## where k in {1,..,d}, i in {1,..,n} [n = length(lx)]
    B <- l.a.k + outer(k-1, lx)
    method <- match.arg(method)
    res <- switch(method,
                  log.poly = {
                      ## stably compute log(colSums(exp(B))) (no overflow)
                      ## Idea:
                      ## (1) let b_k := log(a_{dk}(theta)) + (k-1)*lx and b_{max} := argmax{b_k}.
                      ## (2) \sum_{k=1}^d a_{dk}(theta)\exp((k-1)*lx) = \sum_{k=1}^d \exp(log(a_{dk}(theta))
                      ##     + (k-1)*lx) = \sum_{k=1}^d \exp(b_k) = \exp(b_{max})*\sum_{k=1}^d
                      ##     \exp(b_k-b_{max})
                      ## (3) => log(\sum...) = b_{max} + log(\sum_{k=1}^d \exp(b_k-b_{max}))
                      ## FIXME: faster if scaled with asymptotics
                      max.B <- apply(B, 2, max)
                      max.B + log(colSums(exp(B - rep(max.B, each=d))))
                  },
                  log1p = {
                      ## use log(1 + sum(<smaller>)) = log1p(sum(<smaller>)),
                      ## but we don't expect it to make a difference
                      ## FIXME: faster if scaled with asymptotics
                      im <- apply(B, 2, which.max) # indices (vector) of maxima
                      n <- length(lx) ; d1 <- d-1L
                      max.B <- B[cbind(im, seq_len(n))] # get max(B[,i])_{i=1,..,n} == apply(B, 2, max)
                      B.wo.max <- matrix(B[unlist(lapply(im, function(j) k[-j])) +
                                           d*rep(0:(n-1), each = d1)], d1, n) # matrix B without maxima
                      max.B + log1p(colSums(exp(B.wo.max - rep(max.B, each = d1))))
                  },
                  poly = {
	              ## brute force ansatz
                      log(colSums(exp(B)))
                  },
              {stop(sprintf("unsupported method '%s' in polyJ", method))})
    if(log) res else exp(res)
}

### ==== other numeric utilities ===============================================

##' Properly compute log(x_1 + .. + x_n) for given log(x_1), .., log(x_n)
##'
##' @title Properly compute the logarithm of a sum
##' @param lx matrix or vector of summands (as log(x_1), .., log(x_n))
##' @return log(x_1 + .. + x_n) [for each row of lx] computed via
##'         log(sum(x)) = log(sum(exp(log(x))))
##'         = log(exp(log(x_max))*sum(exp(log(x)-log(x_max))))
##'         = log(x_max) + log(sum(exp(log(x)-log(x_max)))))
##'         = lx.max + log(sum(exp(lx-lx.max)))
##' @author Marius Hofert
lsum <- function(lx, l.off=apply(lx, 1, max)) {
    if(!is.matrix(lx)) lx <- rbind(lx)
    l.off + log(rowSums(exp(lx - l.off))) # FIXME: for a vector lx, names(lx) == "lx" => maybe remove?
}

##' Properly compute log(-+x_1 -+ .. -+ x_n) for given log(|x_1|), .., log(|x_n|)
##' and sign(x_1), .., sign(x_n)
##'
##' @title Properly compute the logarithm of a sum with signed coefficients
##' @param lxabs matrix or vector of summands (as log(|x_1|), .., log(|x_n|))
##' @param signs corresponding matrix or vector of signs (sign(x_1), .., sign(x_n))
##' @return log(x_1 + .. + x_n) [for each row of lx] computed via
##'         log(sum(x)) = log(sum(signs*exp(log(|x|))))
##'         = log(exp(log(|x|_max))*sum(signs*exp(log(|x|)-log(|x|_max))))
##'         = log(|x|_max) + log(sum(signs*exp(log(|x|)-log(|x|_max)))))
##'         = lxabs.max + log(sum(signs*exp(lxabs-lxabs.max)))
##' @author Marius Hofert
lssum <- function(lxabs, signs, l.off=apply(lxabs, 1, max)){
    if(!is.matrix(lxabs)) lxabs <- rbind(lxabs)
    sum. <- rowSums(signs * exp(lxabs - l.off))
    if(any(sum. <= 0)) stop("lssum found non-positive sums")
    l.off + log(sum.) # FIXME: for a vector lxabs, names(lxabs) == "lxabs" => maybe remove?
}

##' Compute Stirling numbers of the 1st kind
##'
##' s(n,k) = (-1)^{n-k} times
##' the number of permutations of 1,2,…,n with exactly k cycles
##'
##' NIST DLMF 26.8 --> http://dlmf.nist.gov/26.8
##'
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
##'
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
##'
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
    method <- match.arg(method)
    switch(method,
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
##'
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
##' ---
##'     {Li}_{1}(z) = -\ln(1-z)
##'     {Li}_{0}(z) = {z \over 1-z}
##'     {Li}_{-1}(z) = {z \over (1-z)^2}
##'     {Li}_{-2}(z) = {z \,(1+z) \over (1-z)^3}
##'     {Li}_{-3}(z) = {z \,(1+4z+z^2) \over (1-z)^4}
##'     {Li}_{-4}(z) = {z \,(1+z) (1+10z+z^2) \over (1-z)^5}.
##' ---
##' Accordingly the polylogarithm reduces to a ratio of polynomials in
##' z, and is therefore a rational function of z, for all nonpositive
##' integer orders. The general case may be expressed as a finite sum:
##' ---
##' {Li}_{-n}(z) = \left(z \,{\partial \over \partial z} \right)^n \frac{z}{1-z}=
##'     = \sum_{k=0}^n k! \,S(n+1,k+1) \left({z \over {1-z}} \right)^{k+1}
##' \ \ (n=0,1,2,\ldots),
##' ---
##' where S(n,k) are the Stirling numbers of the second
##' kind. Equivalent formulae applicable to negative integer orders are
##' (Wood 1992, § 6):
##' ---
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
    stopifnot(length(s) == 1) # for now
    method <- match.arg(method)
    switch(method,
           "sum" = {
               stopifnot((Mz <- Mod(z)) <= 1, Mz < 1 | Re(s) > 1,
                         n.sum > 99, length(n.sum) == 1)
               if(logarithm)
                   log(z)+log(polynEval((1:n.sum)^-s, z))
               else z*polynEval((1:n.sum)^-s, z)
           },
           "negint-s_Stirling" = {
               stopifnot(s == as.integer(s), s <= 1)
               if(s == 1) return(-log1p(-z)) ## -ln(1 -z)
               r <- z/(1 - z)
               ## if(s == 0) return(r)
               n <- -as.integer(s)
               ## k1 <- seq_len(n+1)# == k+1, k = 0...n
               fac.k <- cumprod(c(1, seq_len(n)))
               S.n1.k1 <- Stirling2.all(n+1) ## == Stirling2(n+1, k+1)
               if(logarithm)
                   log(r)+ log(polynEval(fac.k * S.n1.k1, r))
               else r* polynEval(fac.k * S.n1.k1, r)
           },
           stop("unsupported method ", method))
}

### ==== other NON-numerics ====================================================

##' Function which computes the conditional copula function C(v|u) of v given u
##'
##' @title Conditional copula function
##' @param v parameter v
##' @param u parameter u
##' @param family Archimedean family (name or object)
##' @param theta parameter theta
##' @param log if TRUE log(C(v|u)) is returned
##' @author Marius Hofert
##' Note: for some families, this function makes sense for u == 0 or v == 0
##'       since the corresponding limits can be computed; but not for all
cacopula <- function(v, u, family, theta, log = FALSE) {
    stopifnot(length(u) == length(v), 0 < u, u < 1, 0 < v, v < 1)
    cop <- getAcop(family)
    stopifnot(cop@paraConstr(theta))
    res <- cop@psiDabs(rowSums(cop@psiInv(cbind(u, v), theta)), theta, log = TRUE) +
        cop@psiInvD1abs(u, theta, log = TRUE)
    if(log) res else exp(res)
}

##' Function which computes psiDabs via Monte Carlo
##'
##' @title Computing the absolute value of the generator derivatives via Monte Carlo
##' @param t evaluation points
##' @param family Archimedean family (name or object)
##' @param theta parameter value
##' @param degree order of derivative
##' @param n.MC Monte Carlo sample size
##' @param method different methods
##'        direct: direct MC formula
##'        log:    proper log
##' @param log if TRUE the log of psiDabs is returned
##' @author Marius Hofert
psiDabsMC <- function(t, family, theta, degree=1, n.MC, method=c("direct","log"),
                      log=FALSE){
    V0. <- getAcop(family)@V0(n.MC,theta)
    l.V0. <- degree*log(V0.)
    lx <- -V0. %*% t(t) + l.V0.
    method <- match.arg(method)
    switch(method,
           "direct" = {
               ## old code:
               ## summands <- function(t) mean(exp(-V0.*t + l.V0.))
               ## res <- unlist(lapply(t,summands))
	       ## now:
	       res <- colMeans(exp(lx))
               if(log) log(res) else res
           },
           "log" = {
	       lx.max <- apply(lx, 2, max)
               res <- lx.max + log(rowMeans(exp(t(lx) - lx.max)))
               if(log) res else exp(res)
           },
       {stop(sprintf("unsupported method '%s' in psiDabsMC", method))})
}

##' Function for setting the parameter in an acopula
##'
##' @title Settting the parameter in an acopula
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
	stop("theta (=", format(value), ") does not fulfill paraConstr()")
    x
}


##' Construct "paraConstr" function from an "interval"
##'
##' @title Construct "paraConstr" function from an "interval"
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

##' Get one of our "acopula" family objects by name
##'
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

if(getRversion() < "2.12") ## take the version in R >= 2.12.0 (also export!)
    adjustcolor <- function(col, alpha.f = 1, red.f = 1, green.f = 1,
                            blue.f = 1, offset = c(0,0,0,0),
                            transform = diag(c(red.f, green.f, blue.f, alpha.f)))
{
    stopifnot(length(offset) %% 4 == 0,
              !is.null(d <- dim(transform)), d == c(4,4))
    x <- col2rgb(col, alpha=TRUE)/255
    x[] <- pmax(0, pmin(1,
                        transform %*% x +
                        matrix(offset, nrow=4, ncol=ncol(x))))
    rgb(x[1,], x[2,], x[3,], x[4,])
}
