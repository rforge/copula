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

set.seed(1) # set seed

## ==== timing =================================================================

##' Computes user times for the admissible parameter combinations provided by "taus"
##' @param n number of variates to be generated
##' @param family the (nested) Archimedean family to be timed
##' @param taus the sequence of Kendall's tau to be tested
##' @param digits number of digits for the output
##' @param debug print current state of the timing during timing
##' @return a (tau_0 x tau_1)-matrix with first column indicating the user run 
##'         times for V0 and the other cells the run time for V01 corresponding
##'         to two given taus among "taus" based on the generated V0's
##' @author Marius Hofert
timing <- function(n,family,taus,digits = 3,debug = FALSE){

    ## setup
    f <- function(x) format(round(x,digits),nsmall=digits,trim=TRUE,
                            scientific=FALSE) # function to format output
    l <- length(taus)
    time.V0 <- numeric(l)
    time.V01 <- matrix(,nrow=l,ncol=l)
    copFamily <- getAcop(family)
    thetas <- copFamily@tauInv(taus)

    ## timing (based on user time)
    if(l == 1){
	time.V0 <- system.time(V0 <- copFamily@V0(n,thetas[1]))[1] # only generate V0
	time.V0
    }else{ # l >= 2
	for(i in 1:(l-1)){ # run over all theta0
            time.V0[i] <- system.time(V0 <- copFamily@V0(n,thetas[i]))[1] 
            if(debug) cat("V0:  tau_0 = ",f(taus[i]),"; time = ",f(time.V0[i]),
                          "s\n",sep="")
            for(j in (i+1):l){ # run over all theta1
                time.V01[i,j] <- system.time(V01 <- copFamily@V01(V0,thetas[i],
                                                                  thetas[j]))[1] 
                if(debug) cat("V01: tau_0 = ",f(taus[i]),", tau_1 = ",
                              f(taus[j]),"; time = ",f(time.V01[i,j]),"s\n",
                              sep="")
            }
        }
	time.V0[l] <- system.time(V0 <- copFamily@V0(n,thetas[l]))[1] # case i == l: compute V0
	## create result object
	res <- time.V01
	res[,1] <- time.V0 # use run times for V0 in the first column
	row.names(res) <- taus # use taus as row and column headers
	colnames(res) <- taus
	res
    }
}

## ==== results ================================================================

n <- 10000
taus1 <- c(0.05,(1:3)/10) # for AMH: cannot attain more concordance
taus2 <- c((1:9)/10,0.95) # for Clayton: Problem of "LD" for too small alpha
taus3 <- c(0.05,(1:7)/10) # for Frank: theta1 corresponding to tau > 0.8 is so large that 1-exp(-theta1) is 1 which leads to an error in rLog()
taus4 <- c(0.05,(1:9)/10,0.95) # "whole" range

## ==== AMH ====

timing(n,"AMH",taus1)

##       0.05   0.1   0.2   0.3
## 0.05 0.002 0.003 0.004 0.004
## 0.1  0.002    NA 0.003 0.004
## 0.2  0.002    NA    NA 0.004
## 0.3  0.002    NA    NA    NA

## conclusion:
## V0:  uniformly fast
## V01: uniformly fast

## ==== Clayton (todo: tau0 = 0.05) ====

timing(n,"Clayton",taus2)

##        0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9  0.95
## 0.1  0.003 0.056 0.061 0.059 0.059 0.065 0.063 0.060 0.056 0.052
## 0.2  0.003    NA 0.042 0.041 0.038 0.044 0.043 0.044 0.043 0.043
## 0.3  0.002    NA    NA 0.029 0.029 0.029 0.031 0.030 0.026 0.030
## 0.4  0.003    NA    NA    NA 0.021 0.022 0.021 0.019 0.022 0.023
## 0.5  0.003    NA    NA    NA    NA 0.018 0.019 0.018 0.019 0.018
## 0.6  0.003    NA    NA    NA    NA    NA 0.015 0.016 0.015 0.016
## 0.7  0.002    NA    NA    NA    NA    NA    NA 0.013 0.013 0.013
## 0.8  0.003    NA    NA    NA    NA    NA    NA    NA 0.012 0.012
## 0.9  0.002    NA    NA    NA    NA    NA    NA    NA    NA 0.011
## 0.95 0.003    NA    NA    NA    NA    NA    NA    NA    NA    NA

## conclusion:
## V0:  uniformly fast
## V01: faster for larger parameter; caution: does not work for tau0 = 0.05 (due to "LD" having numerical difficulties)

## ==== Frank (todo: theta1 large) ====

timing(n,"Frank",taus3)

##       0.05   0.1   0.2   0.3   0.4   0.5   0.6    0.7
## 0.05 0.001 0.121 0.153 0.114 0.112 0.113 0.115  0.112
## 0.1  0.001    NA 0.119 0.120 0.121 0.119 0.120  0.120
## 0.2  0.001    NA    NA 0.150 0.149 0.150 0.150  0.228
## 0.3  0.002    NA    NA    NA 0.249 0.246 0.244  0.252
## 0.4  0.002    NA    NA    NA    NA 0.578 0.577  0.571
## 0.5  0.002    NA    NA    NA    NA    NA 2.509  2.406
## 0.6  0.002    NA    NA    NA    NA    NA    NA 21.365
## 0.7  0.002    NA    NA    NA    NA    NA    NA     NA

## conclusion:
## V0:  uniformly fast
## V01: for a given theta0: uniformly fast over theta1; but increasing in theta1

## ==== Gumbel ====

timing(n,"Gumbel",taus4)

##       0.05   0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9  0.95
## 0.05 0.008 0.015 0.009 0.009 0.009 0.010 0.009 0.009 0.009 0.009 0.009
## 0.1  0.008    NA 0.008 0.009 0.009 0.009 0.009 0.009 0.009 0.009 0.009
## 0.2  0.007    NA    NA 0.009 0.009 0.009 0.005 0.009 0.009 0.010 0.010
## 0.3  0.008    NA    NA    NA 0.009 0.009 0.009 0.010 0.009 0.009 0.009
## 0.4  0.008    NA    NA    NA    NA 0.009 0.009 0.009 0.008 0.008 0.009
## 0.5  0.005    NA    NA    NA    NA    NA 0.009 0.009 0.009 0.010 0.010
## 0.6  0.008    NA    NA    NA    NA    NA    NA 0.009 0.009 0.009 0.009
## 0.7  0.008    NA    NA    NA    NA    NA    NA    NA 0.009 0.010 0.009
## 0.8  0.008    NA    NA    NA    NA    NA    NA    NA    NA 0.005 0.009
## 0.9  0.008    NA    NA    NA    NA    NA    NA    NA    NA    NA 0.010
## 0.95 0.008    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA

## conclusion:
## V0:  uniformly fast
## V01: uniformly fast

## ==== Joe ====

timing(n,"Joe",taus4)

##       0.05   0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9  0.95
## 0.05 0.001 0.368 0.361 0.409 0.362 0.376 0.368 0.371 0.374 0.369 0.377
## 0.1  0.001    NA 0.365 0.370 0.374 0.377 0.382 0.378 0.388 0.387 0.381
## 0.2  0.002    NA    NA 0.385 0.398 0.411 0.417 0.475 0.440 0.438 0.445
## 0.3  0.004    NA    NA    NA 0.416 0.454 0.483 0.521 0.545 0.551 0.561
## 0.4  0.004    NA    NA    NA    NA 0.560 0.725 0.759 0.878 1.001 0.898
## 0.5  0.005    NA    NA    NA    NA    NA 0.898 1.008 1.201 1.340 1.371
## 0.6  0.005    NA    NA    NA    NA    NA    NA 1.217 1.509 1.883 1.987
## 0.7  0.006    NA    NA    NA    NA    NA    NA    NA 1.686 2.356 2.579
## 0.8  0.006    NA    NA    NA    NA    NA    NA    NA    NA 2.197 2.739
## 0.9  0.006    NA    NA    NA    NA    NA    NA    NA    NA    NA 1.562
## 0.95 0.006    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA

## conclusion:
## V0:  (almost) uniformly fast (at least fast enough)
## V01: increasing in theta0 and theta1
## V01 depends on theta0 and theta1 [but is bounded due to approx-parameter]









