## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
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


### Model selection for copulas based on (k-)cross-validation #####################


##' @title Model selection for copulas based on (k-)cross-validation
##' @param copula object of type 'copula' representing the copula to be evaluated
##'        as model
##'        (if necessary, parameters will be used as starting values for fitCopula())
##' @param x (n, d)-matrix containing the data
##' @param k number of blocks for cross-validation
##'        if NULL, leave-one out cross-validation is performed
##' @param verbose logical indicating whether a progress bar is shown
##' @return "cross validated log-likelihood"
##' @author Ivan Kojadinovic
xvCopula <- function(copula, x, k=NULL, verbose=TRUE, ...)
{
    ## checks
    stopifnot(is(copula, "copula"))
    if(!is.matrix(x))
    {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    stopifnot(is.numeric(x), (d <- ncol(x)) > 1, (n <- nrow(x)) > 0, dim(copula) == d)
    k <- if (is.null(k)) n else as.integer(k)
    stopifnot(k >= 2L, n %/% k >= 1)

    ## setup progress bar
    if(verbose) {
	pb <- txtProgressBar(max=k, style=if(isatty(stdout())) 3 else 1)
	on.exit(close(pb)) # on exit, close progress bar
    }

    ## cross-validation
    xv <- 0
    p <- n %/% k # size of all blocks except maybe the last one
    r <- n - k * p # remaining number of lines
    if (r > 0) k <- k + 1L # to account for last block
    m <- rep(p, k) # sizes of blocks
    if (r > 0) m[k] <- r # size of last block
    b <- c(0, cumsum(m)) # 0 + ending line of each block
    v <- matrix(NA, max(p, r), d) # points where copula density will be evaluated

    ## for each block
    for (i in 1:k)
    {
        sel <- (b[i] + 1):b[i+1] # m[i] lines of current block
        ## estimate copula from all lines except those in sel
        u <- pobs(x.not.s <- x[-sel, , drop=FALSE])
        copula <- fitCopula(copula, u, method = "mpl",
                            estimate.variance=FALSE, ...)@copula
        imi <- seq_len(m[i]) # 1:m[i]
        v.i <- v[imi, , drop=FALSE] # (for efficiency)
        x.sel <- x[sel, , drop=FALSE]
        ## points where copula density will be evaluated
        for (j in 1:d)
            v.i[,j] <- ecdf(x.not.s[,j])(x.sel[,j])
        nmi <- n - m[i] # == nr. of obs. in x.not.s
        ## rescale v.i[,] to *inside* (0, 1) to avoid values 0 or 1
        ## could map  k/n |--> (k+1)/(n+2) --- or rather the same as ppoints():
        a <- if(nmi <= 10) 3/8 else 1/2 ## k/n |--> (k-a)/(n+1-2a) { = (k-1/2)/n for n >= 10}
        v.i <- (v.i * nmi - a) / (nmi + 1-2*a)
        ## cross-validation for block i
        xv <- xv + mean(dCopula(v.i, copula, log = TRUE))

        ## update progress bar
        if(verbose)
            setTxtProgressBar(pb, i)
    }

    xv / k * n
}
