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

#### Timing for the implemented nested Archimedean copulas

## ==== timing =================================================================

##' Computes user times for the admissible parameter combinations provided by "taus"
##' @param n number of variates to be generated
##' @param family the (nested) Archimedean family to be timed
##' @param taus the sequence of Kendall's tau to be tested
##' @param digits number of digits for the output
##' @param verbose print current state of the timing during timing
##' @return a (tau_0 x tau_1)-matrix with first column indicating the user run
##'         times for V0 and the other cells the run time for V01 corresponding
##'         to two given taus among "taus" based on the generated V0's
##' @author Marius Hofert
timing <- function(n,family,taus,digits = 3,verbose = FALSE){

    ## setup
    f <- function(x) format(round(x,digits),nsmall=digits,trim=TRUE,
                            scientific=FALSE) # function to format output
    l <- length(taus)
    time.V0 <- numeric(l)
    time.V01 <- matrix(,nrow=l,ncol=l)
    copFamily <- getAcop(family)
    thetas <- copFamily@tauInv(taus)

    ## timing (based on user time)
    if(l == 1) {
        system.time(V0 <- copFamily@V0(n,thetas[1]))[1] # only generate V0
    } else { # l >= 2
	for(i in seq_along(thetas)) { # run over all theta0
            time.V0[i] <- system.time(V0 <- copFamily@V0(n,thetas[i]))[1]
            if(verbose) cat("V0:  tau_0 = ",f(taus[i]),"; time = ",f(time.V0[i]),
                          "s\n",sep="")
            if(i < l) for(j in (i+1):l) { # run over all theta1
                time.V01[i,j] <-
                    system.time(V01 <- copFamily@V01(V0,thetas[i], thetas[j]))[1]
                if(verbose) cat("V01: tau_0 = ",f(taus[i]),", tau_1 = ",
                              f(taus[j]),"; time = ",f(time.V01[i,j]),"s\n",
                              sep="")
            }
        }
	## create result object
	res <- time.V01
	res[,1] <- time.V0 # use run times for V0 in the first column
	row.names(res) <- taus # use taus as row and column headers
	colnames(res) <- taus
	res
    }
}
