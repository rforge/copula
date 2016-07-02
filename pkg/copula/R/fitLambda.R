## Copyright (C) 2016 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
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



##' @title Non-parametric Estimators of the Matrix of Tail Dependence Coefficients
##' @param u An (n, d)-data matrix of pseudo-observations
##' @param p A probability (in (0,1)) used as 'cut-off' point
##' @param lower.tail A logical indicating whether the lower or upper tail dependence
##'        coefficient is to be computed
##' @param verbose A logical indicating whether a progress bar is displayed
##' @return Estimate of the (matrix of) tail dependence coefficients
##'         (depending on p)
##' @author Marius Hofert
##' @note - See Jaworksi et al. (2009, p. 231) and Schmid and Schmidt (2007)
##'         for the method implemented below. We use this only for d = 2
##'       - If further methods are implemented, introduce 'method' argument
##'         with name "Schmid.Schmidt" for the method below.
fitLambda <- function(u, p, lower.tail = TRUE, verbose = FALSE)
{
    ## Checking
    if(!is.matrix(u)) u <- rbind(u, deparse.level=0L)
    d <- ncol(u)
    stopifnot(0 <= u, u <= 1, 0 <= p, p <= 1,
              length(p) == 1, d >= 2, is.logical(lower.tail))

    ## Compute Lambda
    Lam <- diag(1, nrow = d)
    if(verbose) {
        pb <- txtProgressBar(max = d, style = if(isatty(stdout())) 3 else 1)
        on.exit(close(pb))
    }
    if(lower.tail) {
        M <- matrix(pmax(0, p-u), ncol = d)
        for(i in 1:(d-1)) {
            for(j in (i+1):d)
                Lam[i,j] <- mean(apply(M[, c(i, j)], 1, prod)) # \hat{Lambda}_{ij}
            if(verbose) setTxtProgressBar(pb, i) # update progress bar
        }
        int.over.Pi <- (p^2 / 2)^2
        int.over.M <- p^3 / 3
    } else { # upper tail-dependence coefficient
        M <- matrix(pmin(p, 1-u), ncol = d)
        for(i in 1:(d-1)) {
            for(j in (i+1):d)
                Lam[i,j] <- mean(apply(M[, c(i, j)], 1, prod)) # \hat{Lambda}_{ij}
            if(verbose) setTxtProgressBar(pb, i) # update progress bar
        }
        int.over.Pi <- (p*(2-p)/2)^2
        int.over.M <- (1-(1-p)^2 * (1+2*p))/3
    }
    Lam <- (Lam - int.over.Pi) / (int.over.M - int.over.Pi) # proper scaling
    if(d == 2) {
        min(max(Lam[1,2], 0), 1)
    } else {
        ## Sanity adjustment
        Lam[Lam < 0] <- 0
        Lam[Lam > 1] <- 1
        ## Symmetrize
        Lam. <- Lam + t(Lam)
        diag(Lam.) <- rep.int(1, d)
        Lam.
    }
}

