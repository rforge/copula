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
    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    stopifnot(is.numeric(x), (d <- ncol(x)) > 1, (n <- nrow(x)) > 0, dim(copula) == d)
    k <- if (is.null(k)) n else as.integer(k)
    stopifnot(k >= 1L, n %/% k >= 1)

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
        sel <- (b[i] + 1):b[i+1] # lines of current block
        ## estimate copula from all lines except those in sel
        u <- pobs(x.not.s <- x[-sel, , drop=FALSE])
        copula <- fitCopula(copula, u, estimate.variance=FALSE, ...)@copula
        imi <- seq_len(m[i])# 1:m[i]
        x.sel <- x[sel, , drop=FALSE]
        ## points where copula density will be evaluated
        for (j in 1:d)
            v[imi,j] <- ecdf(x.not.s[,j])(x.sel[,j])
        nmi <- n - m[i]
### FIXME: really rescale *all* v, instead of just v[imi, ]  ????  :
        v <- v * nmi / (nmi + 1L) # rescale to avoid 1
        v[v==0] <- 1 / (nmi + 1L) # to avoid 0
### END(FIXME)
        ## cross-validation for block i
        xv <- xv + mean(dCopula(v[imi, , drop=FALSE], copula, log = TRUE))

        ## update progress bar
        if(verbose)
            setTxtProgressBar(pb, i)
    }

    xv / k * n
}
