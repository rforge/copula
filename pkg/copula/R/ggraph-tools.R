## Copyright (C) 2012 Marius Hofert and Martin Maechler
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


### Tools for graphical gof tests based on pairwise Rosenblatt trafo ###########

### FIXME  replace by cCopula() in the long run

##' Conditional copula function C(u2|u1) of u2 given u1
##'
##' @title Bivariate ("2") Conditional Copula function
##' @param u2 data vector (numeric(n))
##' @param u1 data vector (numeric(n))
##' @param family copula family
##' @param theta parameter (for ACs; for elliptical copulas, its rho)
##' @param ... additional args (e.g., df for t copulas)
##' @return C(u2|u1)
##' @author Marius Hofert and Martin Maechler
##' Note: used in Hofert and Maechler (2013)
ccop2 <- function(u2, u1, family, theta, ...) {
    stopifnot(length(u1)==length(u2))
    switch(copFamilyClass(family),
           "ellipCopula"={
               switch(family,
                      "normal"={
                          pnorm((qnorm(u2)-theta*qnorm(u1))/sqrt(1-theta^2))
                      },
                      "t"={
                          stopifnot(hasArg(df))
                          df <- list(...)$df
                          qt.u1 <- qt(u1, df=df)
                          mu <- theta * qt.u1
                          sigma <- sqrt((df+qt.u1^2)*(1-theta^2)/(df+1))
                          pt((qt(u2, df=df)-mu)/sigma, df=df+1)
                      },
                      stop("not yet supported family"))
           },
           "outer_nacopula"={
               ## cacopula(u, cop=onacopulaL(family, list(theta, 1:2)))
               cop <- getAcop(family)
               u <- cbind(u1, u2)
               psiI <- cop@iPsi(u, theta=theta)
	       exp(cop@absdPsi(rowSums(psiI), theta=theta, log=TRUE) -
		   cop@absdPsi(psiI[,1], theta=theta, log=TRUE))
           },
           stop("family ", family, " not yet supported"))
}

##' Compute pairwise Rosenblatt-transformed variables C(u[,i]|u[,j])) for the given
##' matrix u
##'
##' @title Compute pairwise Rosenblatt-transformed variables
##' @param u (n, d)-matrix (typically pseudo-observations, but also perfectly
##'        simulated data)
##' @param cop copula object used for the Rosenblatt transform (H_0;
##'        either outer_nacopula or ellipCopula)
##' @param ... additional arguments passed to ccop2
##' @return (n,d,d)-array cu.u with cu.u[,i,j] containing C(u[,i]|u[,j]) for i!=j
##'         and u[,i] for i=j
##' @author Marius Hofert
##' Note: used in Hofert and Maechler (2013)
pairwiseCcop <- function(u, cop, ...)
{
    if(!is.matrix(u)) u <- rbind(u)
    stopifnot((d <- ncol(u)) >= 2, 0 <= u, u <= 1,
	      d == dim(cop))

    ## 1) determine copula class and compute auxiliary results
    cls <- copClass(cop)
    family <- copFamily(cop) # determine copula family
    switch(cls,
           "ellipCopula"={
               ## build correlation matrix from vector of copula parameters
               P <- diag(1, nrow=d)
	       rho <- cop@parameters
	       P[lower.tri(P)] <- if(family == "normal" || cop@df.fixed)
		   rho else rho[-length(rho)]
               P <- P + t(P)
               diag(P) <- rep.int(1, d)
           },
           "outer_nacopula"={
               ## build "matrix" of dependence parameters
               P <- nacPairthetas(cop)
           },
           stop("not yet supported copula object"))

    ## 2) compute pairwise C(u_i|u_j)
    n <- nrow(u)
    cu.u <- array(NA_real_, dim=c(n,d,d), dimnames=list(C.ui.uj=1:n, ui=1:d, uj=1:d))
    ## cu.u[,i,j] contains C(u[,i]|u[,j]) for i!=j and u[,i] for i=j
    for(i in 1:d) { # first index C(u[,i]|..)
        for(j in 1:d) { # conditioning index C(..|u[,j])
	    cu.u[,i,j] <- if(i==j) u[,i] else
		ccop2(u[,i], u[,j], family, theta=P[i,j], ...)
        }
    }
    cu.u
}


##' Compute matrix of pairwise tests for independence
##'
##' @title Compute matrix of pairwise tests for independence
##' @param cu.u (n,d,d)-array as returned by \code{pairwiseCcop()}
##' @param N number of simulation \code{N} for \code{\link{indepTestSim}()}
##' @param verbose logical indicating if and how much progress info should be printed.
##' @param ... additional arguments passed to indepTestSim()
##' @return (d,d)-matrix of lists with test results (as returned by indepTest())
##'         for the pairwise tests of independence
##' @author Marius Hofert
##' Note: used in Hofert and Maechler (2013)
pairwiseIndepTest <- function(cu.u, N=256, verbose=TRUE, ...)
{
    ## 1) simulate test statistic under independence
    stopifnot(length(dim. <- dim(cu.u)) == 3L)
    stopifnot(dim.[2]==dim.[3])
    n <- dim.[1]
    d <- dim.[2]
    if(verbose) cat(sprintf("pairwiseIndepTest( (n=%d, d=%d)): indepTestSim(*, N=%d) .. ",
                            n,d, N))
    indepTestDistObj <- indepTestSim(n, p=2, m=2, N=N, ...)
    if(verbose) cat(" *Sim  done\n")

    ## 2) compute matrix of pairwise p-values for test of independence
    p <- matrix(list(), nrow=d, ncol=d)
    for(j in 1:d) { # column
        if(verbose) if(verbose >= 2) cat("j = ", j, ";  i =", sep="") else cat(j, "")
        uj <- cu.u[,,j] # (n x d)
        for(i in 1:d) { # row
            if(verbose >= 2) cat(i,"")
            p[i,j][[1]] <- if(j==i) list(fisher.pvalue = NA)
            else indepTest(cbind(uj[,j], uj[,i]), d=indepTestDistObj)
        }
        if(verbose >= 2) cat("\n")
    }
    p
}

##' Extract p-values from a matrix of indepTest objects
##'
##' @title Extract p-values from a matrix of indepTest objects
##' @param piTest matrix of indepTest objects
##' @return matrix of p-values
##' @author Marius Hofert, Martin Maechler
##' Note: - Kojadinovic, Yan (2010) mention that "fisher.pvalue" performs best;
##'         In d=2 (as we use in pairwiseIndepTest) all three methods are equal.
##'       - used in Hofert and Maechler (2013)
pviTest <- function(piTest){
    matrix(vapply(piTest, function(x) x$fisher.pvalue, numeric(1)),
           ncol=ncol(piTest))
}

##' Computing a global p-value
##'
##' @title Computing a global p-value
##' @param pvalues (matrix of pairwise) p-values
##' @param method vector of methods for p-value adjustments (see ?p.adjust.methods)
##' @param globalFun function determining how to compute a global p-value from a
##'        matrix of pairwise adjusted p-values
##' @return global p-values for each of the specified methods (of how to adjust the
##'         pairwise p-values)
##' @author Marius Hofert
##' Note: used in Hofert and Maechler (2013)
gpviTest <- function(pvalues, method=p.adjust.methods, globalFun=min){
    pvalues <- pvalues[!is.na(pvalues)] # vector
    if(all(c("fdr","BH") %in% method))## p.adjust():  "fdr" is alias for "BH"
	method <- method[method != "fdr"]
    sapply(method, function(meth) globalFun(p.adjust(pvalues, method=meth)))
}

## build global pairwise independent test result string
gpviTest0 <- function(pvalues) {
    pvalues <- pvalues[!is.na(pvalues)] # vector
    c("minimum" = min(pvalues),
      "global (Bonferroni/Holm)" = min(p.adjust(pvalues, method="holm")))
}

## string formatter of gviTest0()
gpviString <- function(pvalues, sep="   ", digits=2) {
    pv <- gpviTest0(pvalues)
    paste0("p-values:", sep,
           paste(paste(names(pv), sapply(pv, format.pval, digits=digits),
                       sep=": "), collapse=paste0(";",sep)))
}


### Tools for graphical gof tests for copulas with radial parts ################

##' Compute pseudo-observations of the radial part and the uniform distribution
##' on the unit sphere (for elliptical copulas) or on the unit simplex (for
##' Archimedean copulas)
##'
##' @title Compute pseudo-observations of the radial part and uniform distribution
##' @param x (n, d)-matrix of data
##' @param do.pobs logical indicating whether pseudo-observations should be computed
##' @param method method slot for different copula classes
##' @param ... additional arguments
##' @return list of two components:
##'         R: n-vector of pseudo-observations of the radial part
##'         U: (n, d)-matrix of pseudo-observations of the uniform distribution on
##'            the unit sphere (for method="ellip") or unit simplex (for method="archm")
##' @author Marius Hofert
##' Note: - the following correct (but unknown!) functions can be provided via '...':
##'         "ellip": qGg, the quantile function of the function G_g in Genest, Hofert, Neslehova (2013)
##'         "archm": iPsi, the inverse of the (assumed) generator
##'       - used in Genest, Hofert, Neslehova (2013)
RUpobs <- function(x, do.pobs=TRUE, method=c("ellip", "archm"), ...)
{
    ## check
    if(!is.matrix(x)) x <- rbind(x, deparse.level=0L)
    method <- match.arg(method)

    u <- if(do.pobs) pobs(x) else x
    switch(method,
           "ellip"={
               ## estimate the standardized dispersion matrix with entries
               ## P_{jk} = \Sigma_{jk}/sqrt{\Sigma_{jj}\Sigma_{kk}}
               ## and compute the inverse of the corresponding Cholesky factor
               P <- nearPD(sin(cor(x, method="kendall")*pi/2), corr=TRUE)$mat
               P <- as.matrix(P) # otherwise, chol() fails
               A. <- chol(P) # upper triangular R such that R^TR = P; note: L = t(A.) s.t. LL^T = P
               A.inv <- solve(A.) # A.^{-1} = (A^{-1})^T

               ## compute Ys
               y <- if(hasArg(qGg)) { # if qGg() has been provided
                   qGg <- list(...)$qGg
                   apply(u, 2, qGg)
               } else { # estimate via empirical quantiles
                   n <- nrow(x)
                   mu <- matrix(rep(colMeans(x), each=n), nrow=n)
                   sd <- matrix(rep(apply(x, 2, sd), each=n), nrow=n)
                   x. <- (x-mu)/sd # standardized data
                   sapply(1:ncol(u), function(j) quantile(x.[,j], probs=u[,j]))
               }

               ## compute Zs and Rs (pseudo-observations of the radial part)
               z <- y %*% A.inv # efficient computation of A^{-1}Y which works with upper triangular matrix A.
               r <- sqrt(rowSums(z^2))

               ## return Rs and Us (pseudo-observations of the uniform distribution on the unit sphere)
               list(R=r, U=z/r)
           },
           "archm"={
               if(hasArg(iPsi)) { # if iPsi() has been provided
                   iPsi <- list(...)$iPsi
                   iPsiu <- iPsi(u)
                   r <- rowSums(iPsiu)
                   list(R=r, U=iPsiu/r)
               } else { # estimate iPsi via inverting the estimator of psi of Genest, Neslehova, Ziegel (2011; Algorithm 1; Section 4.3)
                   ## compute w (slightly adjusted pseudo-observations)
                   n <- nrow(u)
                   W <- vapply(1:n, function(k) sum( colSums(t(U) < u[k,])==d ) / (n+1), NA_real_)
                   ## W's, p's
                   Wtab <- table(W)
                   m <- length(w <- as.numeric(dimnames(Wtab)$x)) # contains w's (= sort(unique(W)))
                   Wtab. <- as.numeric(Wtab)
                   p <- Wtab./n # relative frequency
                   ## solve the system of equations (18) in Genest, Neslehova, Ziegel (2011) to find s
                   if(m < 2) stop("quantity 'm' must be >= 2")
                   s <- numeric(m-1)
                   d <- ncol(u)
                   s[m-1] <- 1-(w[2]/p[m])^(1/(d-1))
                   if(m >= 3){ # note: could be solved explicitly for d==2
                       for(k in (m-2):1){
                           f <- function(x){
                               s. <- c(s[(m-1):(k+1)], x) # at least two elements; the last one is the running variable
                               s.d1 <- (1-cumprod(s.))^(d-1)
                               sum(s.d1*p[(k+1):m])-w[m-k+1]
                           }
                           s[k] <- uniroot(f, interval=c(0,1))$root
                       }
                   }
                   ## construct r's from s's
                   r <- rep(1, m)
                   if(m >= 2) for(k in (m-1):1) r[k] <- r[k+1]*s[k]
                   ## repeat accordingly
                   r <- vapply(1:m, FUN=function(k) rep(r[k], Wtab.[k]), FUN.VALUE=numeric(Wtab.[k]))

                   ## compute
                   iPsihu <- matrix(NA, nrow=n, ncol=d)
                   U <- iPsihu/r # hat{psi}^{-1}(u)/r
                   warning("pseudo-observations 'U' of the uniform distribution on the unit simplex not yet implemented") # TODO (requires inverting the estimator of Genest, Neslehova, Ziegel (2011)) + test this code

                   ## return
                   list(R=r, U=U)
               }
           },
           stop("wrong method"))
}

##' Compute supposedly Beta distributed observations from supposedly uniformly
##' distributed observations on the unit sphere
##'
##' @title Compute supposedly Beta distributed observations from supposedly
##'        uniformly distributed observations on the unit sphere
##' @param u (n, d)-matrix of supposedly uniformly distributed (row) vectors
##'        on the unit sphere in IR^d
##' @return (n, d-1)-matrix where the kth column contains supposedly
##'         Beta(k/2, (d-k)/2)-distributed values
##' @author Marius Hofert
##' Note: - see Li, Fang, Zhu (1997); suggestion: take k~=d/2
##'       - used in Genest, Hofert, Neslehova (2013)
gofBTstat <- function(u){
    if (!is.matrix(u)) u <- rbind(u)
    u2 <- u^2
    t(apply(u2, 1, cumsum))[, -ncol(u)]/rowSums(u2)
}

