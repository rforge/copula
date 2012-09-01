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

user <- Sys.getenv("USER")
switch(user,
       "mhofert"={ # account Marius Hofert (local)
           rm(list=objects()) # remove previously loaded objects
           setwd("~/R/MMMH/parallel/")
       },
       "maechler"={
           setwd("~/R/MM/Pkg-ex/copula")
           Mlibrary(copula)
       },
       stop("unsupported user"))

options(warn=1) # print warnings as they occur

## load packages
stopifnot(require(copula))



### Simulation  a la  Genest et al (2009) -- of only one case 
### -------------------------------------
gofC <- function(copula, H0copula, n, N, method, seed=NULL,
                 verbose=FALSE, ...) {
    if(!is.null(seed)) set.seed(seed)
    T <- system.time(
        pV <- gofCopula(copula, x = rCopula(n, H0copula), N=N,
                        method=method, verbose=verbose, ...)$pvalue)
    c(pvalue = pV, T)# or rather list?
}


require("parallel")
(nC <- detectCores())

## As Genest et al.
n <- 150
dim <- 2
N <- 1000

method <- "SnC"

nRep <- 1000
## Testing:
nRep <- 20

## we set seed *INSIDE* set.seed(17)
(thCl <- iTau(claytonCopula(), 0.5)) # 2
rr <- mclapply(1:nRep, function(i)
	       gofC(normalCopula(dim=dim), n = n, N = N,
		    H0copula = claytonCopula(thCl, dim=dim),
		    seed = i, method = method),
               mc.set.seed=FALSE,
               mc.cores = detectCores())
save(rr, n,dim,N,method,  file = "gof-sim_rr.rda")
