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

##' Find a root
##' @param f function
##' @param interval interval
##' @param lower lower endpoint
##' @param upper upper endpoint
##' @param f.lower function value at lower endpoint
##' @param f.upper function value at upper endpoint
##' @param Sig sign of f.upper
##' @param tol tolerance
##' @param maxiter maximal number of iterations
##' @param trace number determining tracing
##' @author Martin Maechler (from Martin's package "nor1mix")
safeUroot <- function(f, interval, lower = min(interval), upper = max(interval),
		      f.lower = f(lower), f.upper = f(upper),
		      Sig = sign(f.upper),
		      tol = .Machine$double.eps^0.25, maxiter = 1000, trace = 0,
                      ...)
{
    if(trace >= 2)
	cat(sprintf("search in [%g,%g]\n", lower, upper))

    ## make sure we have Sig*f(lower) < 0 and Sig*f(upper) > 0:
    delta.r <- 0.01*max(1e-7, abs(lower))
    f.lo <- f.lower
    while(Sig*f.lo > 0) {
	f.lo <- f(lower <- lower - delta.r)
	if(trace)
	    cat(sprintf(" .. modified lower: %g\n",lower))
	delta.r <- 2 * delta.r
    }
    delta.r <- 0.01*max(1e-7, abs(upper))
    f.up <- f.upper
    while(Sig*f.up < 0) {
	f.up <- f(upper <- upper + delta.r)
	if(trace)
	    cat(sprintf(" .. modified upper: %g\n",upper))
	delta.r <- 2 * delta.r
    }

    ## here, we require R >= 2.6.0 with the improved uniroot():
    uniroot(f, lower=lower, upper=upper,
	    f.lower = f.lo, f.upper = f.up,
	    tol=tol, maxiter=maxiter, ...)
}
