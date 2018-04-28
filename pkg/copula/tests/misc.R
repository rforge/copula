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

### Miscellaneous tests for 'copula' package
### -------------------      ======

require(copula)
source(system.file("Rsource", "utils.R",     package="copula", mustWork=TRUE))
##-> assertError(), assert.EQ(), ... showProc.time()

if(!dev.interactive(orNone=TRUE)) pdf("misc-tests_copula.pdf")


## Testing  chiPlot() and Kplot()  --> ../R/graphics.R
##          -------       -----        ~~~~~~~~~~~~~~~
X <- cbind(c(-2.224, -1.538, -0.807, 0.024, 0.052,  1.324),
           c(0.431,   1.035,  0.586, 1.465, 1.115, -0.847))
## For now,  ":::" because they are not exported yet:
chiP <- copula:::chiPlot(X, main = quote(list(chiPlot(X), X %in% R^{6 %*% 2})))
Kplt <- copula:::Kplot  (X, main = "Kplot(X)")
##
stopifnot(all.equal(Kplt[,"H"], c(0, 0, 2, 2, 6, 6)/10),
          all.equal(Kplt[,"W"], tolerance = 5e-6,
                    c(0.0383659, 0.0924076, 0.163329, 0.255939, 0.381122, 0.568837)),
          all.equal(chiP, tolerance = 2e-6,
                    cbind(H= c(0, 0.2, 0.2, 0.6, 0.6, 0),
                          F= c(0,   0.2, 0.4, 0.6, 0.8, 1),
                          G= c(0.2, 0.6, 0.4, 1,   0.8, 0),
                          chi=c(NaN, 0.408248, 1/6, NaN, -0.25, NaN),
                          lambda=c(1, -0.36, 0.04, 1, 0.36, -1))),
          TRUE)


