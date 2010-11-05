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

#### Implementation of function evaluations and random number generation for
#### nested Archimedean copulas

##' Returns the copula density at u
##' Basically, this is a wrapper based on the slot "dAc" of the provided families
##' which chooses some default values of the parameters
##' warning / FIXME: good default values depend on theta, the dimension, and even
##'                  the evaluation point 
##' @param x nacopula
##' @param u argument of the copula x
##' @param theta copula parameter
##' @param MC TRUE if the derivatives of order > 1 should be evaluated by Monte Carlo
##' @param N approximation parameter for MC = FALSE; sample size for MC = TRUE
##' @param log if TRUE, the log-density is evaluated
##' @author Marius Hofert
dnacopula <- function(x,u,theta,MC = FALSE,N,log = FALSE){
    if(x@childCops != list()){
        stop("currently, only Archimedean copulas are provided")
    }
    cop <- x@copula
    switch(cop@name,
	"AMH" = {cop@dAc(u,theta,MC,N = if(MC) 20000 else 1000,log)},
	"Clayton" = {cop@dAc(u,theta,log)},
	"Frank" = {cop@dAc(u,theta,MC,N = if(MC) 20000 else 1000,log)},
	"Gumbel" = {cop@dAc(u,theta,MC,N = if(MC || ncol(u) > 10) 20000 else 5000,log)},
	"Joe" = {cop@dAc(u,theta,MC,N = if(MC) 20000 else 5000,log)}, 
	{stop("wrong Archimedean family provided")})
}

##' Returns the copula value at a certain vector u
##' @param x nacopula
##' @param u argument of the copula x
##' @return x(u)
##' @author Marius Hofert, Martin Maechler
## FIXME: maybe make this applicable to a matrix of u's?
pnacopula <- function(x,u) {
    stopifnot(is.numeric(u), 0 <= u, u <= 1,
	      length(u) >= dim(x))	# can be larger
    C <- x@copula
    th <- C@theta
    ## use u[j] for the direct components 'comp':
    C@psi(sum(unlist(lapply(u[x@comp], C@psiInv, theta=th)),
	      C@psiInv(unlist(lapply(x@childCops, pnacopula, u = u)),
		       theta=th)),
	  theta=th)
}

##' Compute the probability P[l < U <= u]  where U ~ copula x.
##' @param x outer nested archimedean copula
##' @param l d-vector of lower "integration" limits
##' @param u d-vector of upper "integration" limits
##' @return the probability that a random vector following the given copula
##'  falls in the hypercube with lower and upper corner l and u, respectively.
##' @author Marius Hofert, Martin Maechler
setGeneric("prob", function(x, l, u) standardGeneric("prob"))

setMethod("prob", signature(x ="outer_nacopula"),
          function(x, l,u) {
              d <- dim(x)
              ## TODO: maybe allow  l & u to be  k x d matrices
              ##        --> return vector of probabilities of length k
              stopifnot(is.numeric(l), is.numeric(u),
                        length(u) == d, d == length(l),
                        0 <= l, l <= u, u <= 1)
              if(d > 30)
                  stop("prob() for copula dimensions > 30 are not supported (yet)")
              D <- 2^d
              m <- 0:(D - 1)
              ## digitsBase() from package 'sfsmisc' {slightly simplified} :
              ## Purpose: Use binary representation of 0:N
              ## Author: Martin Maechler, Date:  Wed Dec  4 14:10:27 1991
              II <- matrix(0, nrow = D, ncol = d)
              for (i in d:1L) {
                  II[,i] <- m %% 2L + 1L
                  if (i > 1) m <- m %/% 2L
              }
              ## Sign: the ("u","u",...,"u") case has +1; = c(2,2,...,2)
              Sign <- c(1,-1)[1L + (- rowSums(II)) %% 2]
              U <- array(cbind(l,u)[cbind(c(col(II)), c(II))], dim = dim(II))
              sum(Sign * apply(U, 1, pnacopula, x=x))
          })

##' Returns (n x d)-matrix of random variates
##' @param n number of vectors of random variates to generate
##' @param x outer_nacopula
##' @return matrix of random variates
##' @author Marius Hofert, Martin Maechler
rnacopula <- function(n, x, ...)
{
    Cp <- x@copula			# outer copula
    theta <- Cp@theta			# theta for outer copula
    V0 <- Cp@V0(n,theta)		# generate V0's
    childL <- lapply(x@childCops, rnchild, # <-- start recursion
		     theta0=theta,V0=V0,...)
    dns <- length(x@comp)	 # dimension of the non-sectorial part
    r <- matrix(rexp(n*dns), n, dns) # generate the non-sectorial part
    ## put pieces together
    mat <- Cp@psi(r/V0, theta=theta)	# transform
    mat <- cbind(mat, do.call(cbind,lapply(childL, `[[`, "U")))
    ## get correct sorting order:
    j <- c(x@comp, unlist(lapply(childL, `[[`, "indCol")))
    ## extra checks TODO: comment
    stopifnot(length(j) == ncol(mat))
    m <- mat[,order(j)]			# permute data and return
    ## extra checks TODO: comment
    stopifnot(length(dm <- dim(m)) == 2, dm == dim(mat))
    m
}

##' Returns a list with an (n x d)-matrix of random variates and a vector of
##' indices.
##' @param x nacopula
##' @param n number of vectors of random variates to generate
##' @param theta0 parameter theta0
##' @param V0 vector of V0's
##' @return list(U = matrix(*,n,d), indCol = vector of length d)
##' @author Marius Hofert, Martin Maechler
rnchild <- function(x, theta0, V0,...)
{
    n <- length(V0)
    Cp <- x@copula # inner copula
    ## Consistency checks -- for now {comment later} :
    stopifnot(is(Cp, "acopula"), is.numeric(n), n == as.integer(n),
              is.numeric(V0), length(V0) == n, is.numeric(theta0))
    theta1 <- Cp@theta # theta_1 for inner copula
    ## generate V01's (only for one sector since the
    ## recursion in rnacopula() takes care of all sectors):
    V01 <- Cp@V01(V0, theta0,theta1,...)
    childL <- lapply(x@childCops, rnchild, # <-- recursion
                     theta0=theta1, V0=V01,...)
    dns <- length(x@comp)	# dimension of the non-sectorial part
    r <- matrix(rexp(n*dns), n, dns) # generate the non-sectorial part
    ## put pieces together: first own comp.s, then the children's :
    mat <- Cp@psi(r/V01, theta1) # transform
    mat <- cbind(mat, do.call(cbind, lapply(childL, `[[`, "U")))
    ## get correct sorting order:
    j <- c(x@comp, unlist(lapply(childL, `[[`, "indCol")))
    list(U = mat, indCol = j) # get list and return
}

if(FALSE) { # evaluate the following into your R session if you need debugging:
    trace(rnacopula, browser, exit=browser, signature=signature(x ="outer_nacopula"))

    debug(rnchild)
}

##' Constructor for outer_nacopula
##' @param family either character string (short or longer form of
##'	 copula family name) or an "acopula" family object
##' @param nacStructure a "formula" of the form C(th, comp, list(C(..), C(..)))
##' @return a valid outer_nacopula object
##' @author Martin Maechler
onacopula <- function(family, nacStructure) {
    ## , envir = ... , enclos=parent.frame() or
    ## , envir=environment()
    ##
### FIXME: base this on  onacopulaL() -- replacing nacStructure by nacList
### ----- : (1) replacing  C() with list() *AND* by
###         (2) wrapping the 3rd argument with list(.) if it's not already
### Use a *recursive* function like this one (but add the "list(.)" wrapping:
###    repC <- function(e) { sC <- as.symbol("C"); if(identical(e[[1]], sC)) e[[1]] <- as.symbol("list"); if(length(e) == 4) e[[4]] <- repC(e[[4]]); e}
    nacl <- substitute(nacStructure)
    stopifnot(identical(nacl[[1]], as.symbol("C")))
    COP <- getAcop(family)
    nacl[[1]] <- as.symbol("oC")
### does not work ..>>>>>>>>>>> we should use onacopulaL() inside functions! <<<<
    ## needed, e.g., when onacopula() is called from inside a user function:
    ## for(j in 2:3) if(is.language(nacl[[j]]))
    ##     nacl[[j]] <- eval(nacl[[j]], parent.frame())
    ## pframe <- parent.frame()
    mkC <- function(cClass, a,b,c) {
	if(missing(b) || length(b) == 0) b <- integer()
	if(missing(c) || length(c) == 0) c <- list()
	else if(length(c) == 1 && !is.list(c)) c <- list(c)
### does not work ...
	## else if(!is.list(c)) {
        ##     c <- eval(c, pframe) ; stopifnot(is.list(c))
        ## }
	## if(!is.numeric(a)) {
        ##     a <- eval(a, pframe) ; stopifnot(is.numeric(a), length(a) == 1)
        ## }
	## if(!is.numeric(b)) {
        ##     b <- eval(b, pframe) ; stopifnot(is.numeric(b))
        ## }
	else stopifnot(is.list(c))
	stopifnot(is.numeric(a), length(a) == 1, is.numeric(b))
	if(any(sapply(c, class) != "nacopula"))
	    stop("third entry of 'nacStructure' must be NULL or list of 'C(..)' terms")
	new(cClass, copula = setTheta(COP, a),
	    comp = as.integer(b), childCops = c)
    }
    C <- function(a,b,c) mkC("nacopula", a,b,c)
    oC <- function(a,b,c) mkC("outer_nacopula", a,b,c)
    eval(nacl)
}

##' @title Constructor for outer_nacopula - with list() input and *recursively*
##' @param family
##' @param nacList
##' @return
##' @author Martin Maechler
onacopulaL <- function(family, nacList) {
    COP <- getAcop(family)
    mkC <- function(cClass, abc) {
        stopifnot(is.list(abc), 2 <= (nL <- length(abc)))
        if(nL == 2) abc[[3]] <- list() else
        if(nL > 3) stop("nacLists must be of length 3 (or 2)")
	new(cClass, copula = setTheta(COP, abc[[1]]),
	    comp = as.integer(abc[[2]]),
            childCops = lapply(abc[[3]], function(.) mkC("nacopula", .)))
    }
    mkC("outer_nacopula", nacList)
}

