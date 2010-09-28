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

#### Estimation for (nested) Archimedean copulas

##' Computes the maximum likelihood estimator for an Archimedean copula
##' @param u matrix of realizations following the copula
##' @param cop acopula to be estimated
##' @return MLE
##' @author Marius Hofert
mleacopula <- function(u,cop){
    ## remark: this is not optimal yet due to two reasons:
    ## (1) -log(density()) is only optained approximately (which is a problem 
    ## 	   particularly for Gumbel and Joe)
    ## (2) -log(density()) evaluates the setup steps "not depending on theta" 
    ## 	   for every call---this is inefficient
    interval <- switch(cop@name,
                       AMH = {c(0.01,0.99)}, # upper bound corresponds to tau = 0.33
                       Clayton = {c(0.01,4.666667)}, # upper bound corresponds to tau = 0.7
                       Frank = {c(0.01,11.41157)}, # upper bound corresponds to tau = 0.7
                       Gumbel = {c(1.01,3.333333)}, # upper bound corresponds to tau = 0.7
                       Joe = {c(1.01,5.463749)}, # upper bound corresponds to tau = 0.7
                   {stop("wrong family")})
    optimize(cop@mLogDensity,interval,tol=0.001)$minimum
}	

