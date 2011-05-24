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

require(nacopula)
sessionInfo() # will change too often.. but if we need the info

### ==== A faster, more testing version of  demo(estimation.gof)  ==============
##					i.e., ../demo/estimation.gof.R :
##					      ~~~~~~~~~~~~~~~~~~~~~~~~

source(system.file("Rsource", "estim-gof-fn.R", package="nacopula"))
## --> estimation.gof() etc

## Use all available estimation and GOF methods:
(estMeth <- eval(formals(enacopula)$method))
(gofMeth <- eval(formals(gnacopula)$method))

set.seed(1) # set seed

n <- 100 # sample size -- small here for CPU reasons
d <- 5 # dimension
tau <- 0.2 # Kendall's tau

## ==== apply all procedures (to data from AMH) ================================

simFamily <- "AMH"
cop <- getAcop(simFamily)
theta <- cop@tauInv(tau) # true parameter

## start the loop
cat("\n## ==== data from ",simFamily," (n = ",n,", d = ",d,", theta = ",
    format(theta),", tau = ", format(tau),") ====\n\n",sep="")

if(getRversion() <= "2.13")
    source(system.file("Rsource", "fixup-sapply.R", package="nacopula"))

RR <- sapply(estMeth, simplify="array", function(e)
         {
             sapply(gofMeth, simplify="array", function(g)
                    estimation.gof(n, d, simFamily, tau = tau, n.MC = 0,
                                   esti.method = e, gof.method = g))
         })

str(RR)
## Now print RR smartly (well, "to be improved"):
options(digits = 5)

## *Not* the times here:
RR[,c(1:2,4:5),,]

### Now do use the parametric bootstrap for better P-value:
set.seed(11)
##
n <- 64    # small sample size
d <- 5     # dimension
tau <- 0.8 # Kendall's tau
(theta <- copGumbel@tauInv(tau)) # == 5  [true parameter]

## now [2011-04-29] should work:  emle() *did* fail for Frank here:
R2 <- sapply(gofMeth, simplify="array", function(g)
	     estimation.gof(n, d, copGumbel, tau = tau, n.MC = 0,
			    n.bootstrap = 256,
			    esti.method = "mle", gof.method = g,
                            checkFamilies = nacopula:::c_longNames[c("C","F","J", "G")]))
R2
## }

cat('Time elapsed: ', proc.time(),'\n') # for ''statistical reasons''

### Make sure the log-Likelihood demos run: ==================================

demo("logL-vis", package = "nacopula")

cat('Time elapsed: ', proc.time(),'\n') # for ''statistical reasons''
