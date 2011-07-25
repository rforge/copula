## Copyright (C) 2010--2011 Marius Hofert and Martin Maechler
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
##' numerically accurately, for both limits  th -> 0  &  th -> 1
##' @author Martin Maechler
tauAMH <- function(th) {
    if(length(th) == 0) return(numeric(0)) # to work with NULL
    r <- th
    na <- is.na(th)
    lrg <- (th > 0.01) & !na
    f. <- function(t) {
	## 1 - 2*((1-t)*(1-t)*log1p(-t) + t) / (3*t*t)
	r <- t
	r[i1 <- (1-t) == 0] <- 1/3
        t <- t[s <- !i1]
	r[s] <- 1 - 2*((1-t)^2 *log1p(-t) + t) / (3*t*t)
	r
    }
    r[lrg] <- f.(th[lrg])
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
m.opt.retst <- function(V0) {
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
retstablerej <- function(m,V0,alpha) {
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
retstableR <- function(alpha, V0, h=1) {
    n <- length(V0)
    stopifnot(is.numeric(alpha), length(alpha) == 1,
	      0 <= alpha, alpha <= 1) ## <- alpha > 1 ==> cos(pi/2 *alpha) < 0
    ## case alpha==1
    if(alpha == 1 || n == 0) { # alpha == 1 => St corresponds to a point mass at V0 with
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
    stopifnot(is.numeric(alpha), length(alpha) == 1,
	      0 < alpha, alpha <= 1,
	      is.numeric(h), length(h) == 1, h > 0)
    if(alpha == 1 || n == 0) {
	## alpha == 1 => St corresponds to a point mass at V0 with
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
##' @param alpha parameter (1/theta) in (0,1]
##' @param method a character string, one of
##'    "sort":          compute coefficients via exp(log()) pulling out the maximum, and sort
##'    "horner":        uses polynEval()
##'    "direct":        brute force approach
##'    "dsumSibuya":    uses dsumSibuya() - FIXME? allow to specify *which* dsumS.. method
##' @param log boolean which determines if the logarithm is returned
##' @param verbose logical for method == sort
##' @return a_{dk}(theta) = (-1)^{d-k}\sum_{j=k}^d alpha^j * s(d,j) * S(j,k)
##' @author Marius Hofert und Martin Maechler
##' note: this function is known to cause numerical problems, e.g., for d=100, alpha=0.8
coeffG <- function(d, alpha, method = c("sort", "horner", "direct", "dsumSibuya"),
		   log = FALSE, verbose = FALSE)
{
    stopifnot(is.numeric(d), length(d) == 1, d >= 1, length(alpha) == 1,
              0 < alpha, alpha <= 1)
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
	   "dsumSibuya" = {
	       ## coefficients via dsumSibuya
	       ## a_{dk}(theta) = d!/k! * dsumSibuya(d, k, alpha)
	       k <- 1:d
	       ck <- ## c_k := d!/k!
		   if(log) c(0,cumsum(log(d:2)))[d:1]
		   else c(1,cumprod(d:2))[d:1]
	       p <- dsumSibuya(d, k, alpha, log=log)
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
	   stop(sprintf("unsupported method '%s' in coeffG", method))
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
##' @param alpha parameter (1/theta) in (0,1]
##' @param d number of summands, >= 1
##' @param method a string, one of
##'   "default"         uses a combination of the other methods
##'   "pois.direct"     uses ppois directly
##'   "pois"            uses ppois with pulling out max
##'   "stirling"        uses the representation via Stirling numbers and once horner
##'   "stirling.horner" uses the representation via Stirling numbers and twice horner
##'   "sort"            compute the coefficients via exp(log()), pulling out the max and sort
##'   "horner"          uses polynEval
##'   "direct"          brute force approach
##'   "dsumSibuya"            uses dsumSibuya()
##' @param log boolean which determines if the logarithm is returned
##' @return \sum_{k=1}^d  a_{dk}(\theta)  x ^ k
##'       = \sum_{k=1}^d  a_{dk} *     exp(lx*k)
##'  where a_{dk}(theta)
##'       = (-1)^{d-k}\sum_{j=k}^d \theta^{-j} s(d,j) S(j,k)
##'       = (d!/k!)\sum_{l=1}^k (-1)^{d-l} \binom{k}{l}\binom{\alpha l}{d}
##' @author Marius Hofert
polyG <- function(lx, alpha, d, method=c("default", "pois", "pois.direct",
                                "stirling", "stirling.horner", "sort", "horner",
                                "direct", "dsumSibuya"), log=FALSE)
{
    k <- 1:d
    stopifnot(length(alpha)==1, 0 < alpha, alpha <= 1)
    method <- match.arg(method)
    switch(method,
	   "default" =
       {
	   method <- if(d <= 100)
               if(alpha <= 0.54) "stirling" else if(alpha <= 0.77)
                   "pois.direct" else "dsumSibuya"
           else "pois" # slower but more stable, e.g., for d=150
	   polyG(lx=lx, alpha=alpha, d=d, method=method, log=log)
       },
           "pois" =
       {
           ## determine signs of the falling factorials
           signs <- (-1)^k* (2*(floor(alpha*k) %% 2) - 1)

           ## build list of b's
           n <- length(lx)
           x <- exp(lx) ## e^lx = x
           lppois <- outer(d-k, x, FUN=ppois, log.p=TRUE) # a (d x n)-matrix; log(ppois(d-k, x))
           llx <- k %*% t(lx) # also a (d x n)-matrix; k*lx
           labsPoch <- vapply(k, function(j) sum(log(abs(alpha*j-(k-1L)))), NA_real_) # log|(alpha*k)_d|
           lfac <- lfactorial(k)
           ## build matrix of exponents
           lxabs <- llx + lppois + rep(labsPoch - lfac, n) + rep(x, each = d)
	   res <- lssum(lxabs, signs, strict=FALSE)
           if(log) res else exp(res)
       },
           "pois.direct" =
       {
           ## determine signs of (-1)^(d-k)*(alpha*k)_d
           signs <- (-1)^k * (2*(floor(alpha*k) %% 2) - 1)

           ## build coefficients
           xfree <- lchoose(alpha*k,d) + lfactorial(d) - lfactorial(k)
           x <- exp(lx)
           lppois <- outer(d-k, x, FUN=ppois, log.p=TRUE) # (length(x),d)-matrix
           klx <- lx %*% t(k)
           exponents <- exp(t(x+klx)+lppois+xfree) # (d,length(x))-matrix
           res <- as.vector(signs %*% exponents)
           if(log) log(res) else res
       },
           "stirling" =
       {
	   ## implementation of \sum_{k=1}^d a_{dk}(\theta) x^k
	   ## = (-1)^{d-1} * x * \sum_{k=1}^d alpha^k * s(d,k) * \sum_{j=1}^k S(k,j) * (-x)^{j-1}
	   ## = (-1)^{d-1} * x * \sum_{k=1}^d alpha^k * s(d,k) * polynEval(...)
	   ## inner function is evaluated via polynEval
	   x <- exp(lx)
	   s <- Stirling1.all(d) # s(d,1), ..., s(d,d)
	   S <- lapply(k, Stirling2.all) # S[[l]][n] contains S(l,n), n = 1,...,l
	   lst <- lapply(k, function(k.) (-1)^(d-1)*x*alpha^k.*s[k.]*polynEval(S[[k.]],-x))
	   res <- rowSums(matrix(unlist(lst), nrow=length(x)))
	   if(log) log(res) else res
       },
           "stirling.horner" =
       {
	   ## implementation of \sum_{k=1}^d a_{dk}(\theta) x^k
	   ## = (-1)^{d-1} * x * \sum_{k=1}^d alpha^k * s(d,k) * \sum_{j=1}^k S(k,j) * (-x)^{j-1}
	   ## = (-1)^{d-1} * x * \sum_{k=1}^d alpha^k * s(d,k) * polynEval(...)
	   ## polynEval is used twice
	   x <- exp(lx)
	   s <- Stirling1.all(d) # s(d,1), ..., s(d,d)
	   S <- lapply(k, Stirling2.all) # S[[l]][n] contains S(l,n), n = 1,...,l
	   len <- length(x)
           poly <- matrix(unlist(lapply(k, function(k.) polynEval(S[[k.]],-x))), nrow=len) # (len,d)-matrix
           res <- (-1)^(d-1)*alpha*x*unlist(lapply(1:len, function(x) polynEval(s*poly[x,], alpha)))
           if(log) log(res) else res
           ## the following code was *not* faster
           ## poly <- t(sapply(k, function(k.) polynEval(S[[k.]],-x))) # (d,len(x))-matrix
           ## coeff <- if(length(x)==1) t(s*poly) else s*poly
           ## res <- (-1)^(d-1)*alpha*x*apply(coeff, 2, polynEval, x=alpha)
       },
           "sort" =, "horner" =, "direct" =, "dsumSibuya" =
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
           logx <- l.a.dk + k %*% t(lx)
           if(log) {
               ## compute log(colSums(exp(B))) stably (no overflow) with the idea of
               ## pulling out the maxima
               lsum(logx)
           }else colSums(exp(logx))
       },
	   stop(sprintf("unsupported method '%s' in polyG", method))
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
##' @param alpha parameter in (0,1]
##' @return vector of random variates V
##' @author Marius Hofert, Martin Maechler
rSibuyaR <- function(n, alpha) {
    stopifnot((n <- as.integer(n)) >= 0, length(alpha)==1, 0 < alpha, alpha <= 1)
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
##' @param n sample size (has to be numeric, >= 0)
##' @param alpha parameter in (0,1]
##' @return vector of random variates V
##' @author Martin Maechler
rSibuya <- function(n, alpha) {
    .Call(rSibuya_vec_c, n, alpha)
}

##' Probability mass function of a Sibuya(alpha) distribution
##'
##' @title Probability mass function of a Sibuya(alpha) distribution
##' @param x evaluation point [integer]
##' @param alpha parameter alpha
##' @param log boolean which determines if the logarithm is returned
##' @return p_x = choose(alpha, x) * (-1)^(x-1)
##' @author Marius Hofert and Martin Maechler
dSibuya <- function(x, alpha, log=FALSE)
    if(log) lchoose(alpha, x) else abs(choose(alpha, x))

##' Distribution function of a Sibuya(alpha) distribution
##'
##' @title Distribution function of a Sibuya(alpha) distribution
##' @param x evaluation point [integer]
##' @param alpha parameter alpha
##' @param lower.tail if TRUE, probabilities are P[X ≤ x], otherwise, P[X > x]
##' @param log.p boolean which determines if the logarithm is returned
##' @return F(x) = 1 - (-1)^x * choose(alpha-1, x)
##' @author Marius Hofert and Martin Maechler
pSibuya <- function(x, alpha, lower.tail=TRUE, log.p=FALSE)
{
    ## F(x) = 1 - 1/(x*Beta(x,1-alpha)) = 1 - (x*beta(x, 1-alpha))^(-1)
    if(log.p) {
        if(lower.tail) # log(1 - 1/(x*beta(x, 1-alpha)))
            log1p(-1/(x*beta(x, 1-alpha)))
        else ## log(1/(x*beta(x, 1-alpha))) = - log(x * beta(..)) =
            -log(x) - lbeta(x, 1-alpha)
    } else { ## no log
        xb <- 1/(x*beta(x, 1-alpha))
        if(lower.tail) 1 - xb else xb
    }
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

##' Inner probability mass function for a nested Joe copula, i.e. a Sibuya sum
##'
##' @title Inner probability mass function for a nested Joe copula
##' @param x vector (or number) of natural numbers >= n
##' @param n vector (or number) of natural numbers
##' @param alpha parameter in (0,1]
##' @param method method applied
##'        log:      proper log computation based on lssum
##'        direct:   brute-force evaluation of the sum and its log
##'        Rmpfr:    multi-precision
##'        diff:     via forward differences
##'        exp.log:  similar to method = "log", but without *proper/intelligent* log
##' @param log boolean which determines if the logarithm is returned
##' @return p_{xn} = \sum_{j=1}^n choose(n,j)*choose(alpha*j,x)*(-1)^(x-j)
##'         which is a probability mass function in x on IN with generating function
##'         g(z) = (1-(1-z)^alpha)^n
##' @author Marius Hofert and Martin Maechler
##' note: - p_{xn} = 0 for x < n; p_{nn} = alpha^n
##'       - numerically challenging, e.g., dsumSibuya(100, 96, 0.01) < 0 for all methods
dsumSibuya <- function(x, n, alpha,
                       method=c("log", "direct", "Rmpfr", "diff", "exp.log"), log=FALSE)
{
    stopifnot(x == round(x), n == round(n), n >= 1, length(alpha) == 1,
              0 < alpha, alpha <= 1)
    if((l.x <- length(x)) * (l.n <- length(n)) == 0)
	return(numeric())
    if((len <- l.x) != l.n) { ## do recycle to common length
	len <- max(l.x, l.n)
	if(l.x < len)
	    x <- rep(x, length.out = len)
	else ## if(l.n < len)
	    n <- rep(n, length.out = len)
    }
    if(alpha == 1)
	return(x == n)
    ii <- seq_len(len)
    method <- if(missing(method) && is(alpha, "mpfr"))
        "Rmpfr" else match.arg(method)
    switch(method,
	   "log" =
       {
	   ## computes *proper* log based on lssum

	   ## determine the matrix of signs of choose(alpha*j,x)*(-1)^(x-j),
	   ## j in {1,..,m} -- which notably do *not* depend on x !
	   m <- max(n)
	   signs <- unlist(lapply(1:m, function(j) {
	       z <- alpha*j
	       if(z == floor(z)) 0 else (-1)^(j-ceiling(z))
	   }))
	   ## for one pair of x and n:
	   f.one <- function(x,n) {
	       if(x < n) return(-Inf)	# = log(0)
	       j <- seq_len(n)
	       lxabs <- lchoose(n, j) + lchoose(alpha*j, x)
	       ## *NON*-strict -- otherwise need try() :
	       lssum(as.matrix(lxabs), signs[j], strict=FALSE)
	   }
	   S. <- sapply(ii, function(i) f.one(x[i], n[i]))
	   if(log) S. else exp(S.)
       },
	   "direct" =
       {
	   ## brute force evaluation of the sum and its log
	   f.one <- function(x,n) {
	       if(x < n) return(0)
	       j <- seq_len(n)
	       sum(choose(n,j)*choose(alpha*j,x)*(-1)^(x-j))
	   }
	   S <- sapply(ii, function(i) f.one(x[i], n[i]))
	   if(log) log(S) else S
       },
	   "Rmpfr" =
       {
           ## as "direct" but using high-precision arithmetic, where
           ## the precision should be set via alpha = mpfr(*, precBits= .)
	   stopifnot(require(Rmpfr))

           ## FIXME: for one n and many x -- should be made *much* faster!

           if(!is(alpha, "mpfr"))
               alpha <- mpfr(alpha, precB = max(100, min(x, 10000)))
	   mpfr.0 <- mpfr(0, precBits = getPrec(alpha))
	   f.one <- function(x,n) {
	       if(x < n) return(mpfr.0)
	       j <- seq_len(n)
	       sum(chooseMpfr.all(n)*chooseMpfr(alpha*j,x)*(-1)^(x-j))
	   }
	   S <- new("mpfr", unlist(lapply(ii, function(i)
                                          f.one(x[i], n[i]))))
	   as.numeric(if(log) log(S) else S)
       },
           "diff" =
       {
           diff(choose(n:0*alpha, x), differences=n) * (-1)^x
       },
	   "exp.log" =
       {
	   ## similar to method = "log", but without *proper/intelligent* log
	   ## and inefficient due to the signs (old version)
	   f.one <- function(x,n) {
	       if(x < n) return(0)
	       j <- seq_len(n) ## indices of the summands
	       signs <- (-1)^(j+x)
	       ## determine the signs of choose(j*alpha,x) for each component of j
	       to.subtract <- 0:(x-1)
	       sig.choose <-
		   unlist(lapply(j, function(l)
				 prod(sign(l*alpha-to.subtract)) ))
	       signs <- signs*sig.choose
	       binom.coeffs <- exp(lchoose(n,j) + lchoose(j*alpha,x))
	       sum(signs*binom.coeffs)
	   }
	   S <- sapply(ii, function(i) f.one(x[i], n[i]))
	   if(log) log(S) else S
       },
	   ## otherwise
	   stop(sprintf("unsupported method '%s' in dsumSibuya", method)))
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
##' @param alpha parameter alpha ( := 1/theta ) in (0,1]
##' @param d number of summands
##' @param method different methods, can be
##'        "log.poly" intelligent log version
##'        "log1p"    additonally uses log1p
##'        "poly"     brute force log version
##' @param log boolean which determines if the logarithm is returned
##' @return \sum_{k=1}^d a_{d,k}(theta) exp((k-1)*lx) = \sum_{k=0}^{d-1} a_{d,k+1}(theta) x^k
##'         where a_{d,k}(theta) = S(d,k)*(k-1-alpha)_{k-1} = S(d,k)*Gamma((1:d)-alpha)/Gamma(1-alpha)
##' @author Marius Hofert and Martin Maechler
polyJ <- function(lx, alpha, d, method=c("log.poly","log1p","poly"), log=FALSE) {
    stopifnot(length(alpha)==1, 0 < alpha, alpha <= 1)
    ## compute the log of the coefficients a_{dk}(theta)
    if(d > 220) stop("d > 220 not yet supported")# would need Stirling2.all(d, log=TRUE)
    k <- 1:d
    l.a.k <- log(Stirling2.all(d)) + lgamma(k-alpha) - lgamma(1-alpha) # log(a_{dk}(theta)), k = 1,..,d
    ## evaluate polynomial via exp( log(<poly>) )
    ## for this, create a matrix B with (k,i)-th entry B[k,i] = log(a_{dk}(theta)) + (k-1) * lx[i],
    ## where k in {1,..,d}, i in {1,..,n} [n = length(lx)]
    B <- l.a.k + (k-1) %*% t(lx)
    method <- match.arg(method)
    switch(method,
           "log.poly" = {
               ## stably compute log(colSums(exp(B))) (no overflow)
               ## Idea:
               ## (1) let b_k := log(a_{dk}(theta)) + (k-1)*lx and b_{max} := argmax{b_k}.
               ## (2) \sum_{k=1}^d a_{dk}(theta)\exp((k-1)*lx) = \sum_{k=1}^d \exp(log(a_{dk}(theta))
               ##     + (k-1)*lx) = \sum_{k=1}^d \exp(b_k) = \exp(b_{max})*\sum_{k=1}^d
               ##     \exp(b_k-b_{max})
               ## (3) => log(\sum...) = b_{max} + log(\sum_{k=1}^d \exp(b_k-b_{max}))
               if(log) lsum(B) else exp(lsum(B))
           },
           "log1p" = {
               ## use log(1 + sum(<smaller>)) = log1p(sum(<smaller>)),
               ## but we don't expect it to make a difference
               im <- apply(B, 2, which.max) # indices (vector) of maxima
               n <- length(lx) ; d1 <- d-1L
               max.B <- B[cbind(im, seq_len(n))] # get max(B[,i])_{i=1,..,n} == apply(B, 2, max)
               B.wo.max <- matrix(B[unlist(lapply(im, function(j) k[-j])) +
                                    d*rep(0:(n-1), each = d1)], d1, n) # matrix B without maxima
               res <- max.B + log1p(colSums(exp(B.wo.max - rep(max.B, each = d1))))
               if(log) res else exp(res)
           },
           "poly" = {
               ## brute force ansatz
               res <- colSums(exp(B))
               if(log) log(res) else res
           },
       stop(sprintf("unsupported method '%s' in polyJ", method)))
}

##' Circular/Rational function	(1 - x^d)/(1 - x) for x ~~ 1, i.e.,
##' compute (1 - x^d)/(1 - x) = (1 - (1-e)^d) / e   for	 e = 1-x (<< 1) and integer d
##'
##' @title Circular/Rational function  (1 - (1-e)^d) / e  {incl. limit e -> 0}
##' @param e numeric vector in [0, 1]
##' @param d integer (scalar), >= 1
##' @return (1 - (1-e)^d) / e
##' @author Martin Maechler, Date: 25 Jul 2011
circRat <- function(e, d)
{
### TODO (?):  improve "log=TRUE", for e ~= 1:	log(circRat(e, d)) = log1p(-x^d) - log(e)
    stopifnot(length(d) == 1, d == as.integer(d), d >= 1)
    if(d <= 6)
	switch(d,
	       1-0*e, ## <<- d = 1
	       2-e,   ## <<- d = 2: 1 2 1
	       3-e*(3-e),#   d = 3: 1 3 3 1
	       4-e*(6-e*(4-e)),		       ## d = 4: 1 4 6 4 1
	       5-e*(10-e*(10-e*(5-e))),	       ## d = 5: 1 5 10 10 5 1
	       6-e*(15-e*(20-e*(15-e*(6-e))))  ## d = 6: 1 6 15 20 15 6 1
	       )
    else { ## d >= 7 ---------------
	r <- e
	eps <- .Machine$double.eps
	if(any(l1 <- ((d1 <- (d-1)/2)*e < eps)))
	    r[l1] <- d
	if(any(l2 <- !l1 & ((d2 <- (d-2)/3)*(e2 <- e*e) < eps)))
	    r[l2] <- d*(1 - d1*e[l2])
	if(any(l3 <- !l1 & !l2 & ((d3 <- (d-3)/4)*e*e2 < eps)))
	    r[l3] <- d*(1 - d1*e[l3]*(1 - d2*e[l3]))
	## and for the remaining ones, we afford a little precision loss:
	if(any(lrg <- !l1 & !l2 & !l3)) {
	    e <- e[lrg]
	    r[lrg] <- (1 - (1-e)^d)/e
	}
	r
    }
}


### ==== other NON-numerics ====================================================

##' Conditional copula function C(u[,d]|u[,1],...,u[,d-1])
##'
##' @title Conditional copula function
##' @param u (n x d)-matrix of evaluation points (first d-1 colums are conditioned on)
##' @param cop an outer_nacopula
##' @param n.MC Monte Carlo sample size
##' @param log if TRUE the logarithm of the conditional copula is returned
##' @author Marius Hofert
cacopula <- function(u, cop, n.MC=0, log=FALSE) {
    stopifnot(is(cop, "outer_nacopula"))
    if(length(cop@childCops))
	stop("currently, only Archimedean copulas are provided")
    if(!is.matrix(u)) u <- rbind(u)
    stopifnot(0 <= u, u <= 1)
    ## note: for some (but not all) families, cacopula also makes sense on
    ##       the boundaries since the corresponding limits exist
    th <- cop@copula@theta
    stopifnot(cop@copula@paraConstr(th))
    dim. <- dim(u)
    n <- dim.[1]
    d <- dim.[2]
    psiI <- cop@copula@psiInv(u, theta=th)
    arg.denom <- rowSums(psiI[,1:(d-1), drop=FALSE])
    arg.num <- arg.denom + psiI[,d]
    logD <- cop@copula@psiDabs(c(arg.num, arg.denom), theta=th, degree=d-1,
                               n.MC=n.MC, log=TRUE)
    res <- logD[1:n]-logD[(n+1):(2*n)]
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
##'        log:         proper log using lsum
##'        direct:      direct evaluation of the sum
##'        pois.direct: directly uses the Poisson density
##'        pois:        intelligently uses the Poisson density with lsum
##' @param log if TRUE the log of psiDabs is returned
##' @author Marius Hofert
##' Note: psiDabsMC(0) is always finite, although, theoretically, psiDabs(0) may
##'       be Inf (e.g., for Gumbel and Joe)
psiDabsMC <- function(t, family, theta, degree=1, n.MC,
                      method=c("log", "direct", "pois.direct", "pois"),
                      log = FALSE)# not yet: is.log.t = FALSE)
{
    res <- numeric(length(t))
    V <- getAcop(family)@V0(n.MC, theta)
    method <- match.arg(method)
    switch(method,
	   ## the following is not faster than "log":
	   ## "default" = { # basically, use "direct" if numerically not critical and "log" otherwise
	   ##                lx <- -V %*% t(t) + degree*log(V)
	   ##                explx <- exp(lx) # (n.MC, n)-matrix containing the summands
	   ##                explx0 <- explx==0 # can be TRUE due to t == Inf or t finite but too large
	   ##                t.too.large <- unlist(lapply(1:n, function(x) any(explx0))) # boolean vector of length n indicating which column of explx contains zeros
	   ##                r1 <- colMeans(explx[,!t.too.large, drop=FALSE])
	   ##                res[!t.too.large] <- if(log) log(r1) else r1
	   ##                r2 <- lsum(lx[,t.too.large, drop=FALSE] - log(n.MC))
	   ##                res[t.too.large] <- if(log) r2 else exp(r2)
	   ##                res[is.infinite(t)] <- if(log) -Inf else 0
	   ##                res
	   ##            },
           "log" = { # intelligent log
               iInf <- is.infinite(t)
               res[iInf] <- -Inf # log(0)
               if(any(!iInf))
                   res[!iInf] <- lsum(-V %*% t(t[!iInf]) + degree*log(V) - log(n.MC))
               if(log) res else exp(res)
           },
           "direct" = { # direct method
               lx <- -V %*% t(t) + degree*log(V)
               res <- colMeans(exp(lx)) # can be all zeros if lx is too small [e.g., if t is too large]
               if(log) log(res) else res
           },
           "pois.direct" = {
               poi <- dpois(degree, lambda=V %*% t(t))
               res <- factorial(degree)/t^degree * colMeans(poi)
               if(log) log(res) else res
           },
           "pois" = {
               iInf <- is.infinite(t)
               res[iInf] <- -Inf # log(0)
               if(any(!iInf)) {
                   t. <- t[!iInf]
                   lpoi <- dpois(degree, lambda=V %*% t(t.), log=TRUE) # (n.MC, length(t.))-matrix
                   b <- -log(n.MC) + lfactorial(degree) - degree*rep(log(t.), each=n.MC) + lpoi # (n.MC, length(t.))-matrix
                   res[!iInf] <- lsum(b)
               }
               if(log) res else exp(res)
           },
	   stop(sprintf("unsupported method '%s' in psiDabsMC", method)))
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
mkParaConstr <- function(int) {
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
