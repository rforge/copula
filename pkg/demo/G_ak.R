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

library(nacopula)

#### Testing /  exploring  coeffG(), the coefficients  a_k  for
#### the Gumbel copula's generator derivatives and copula density

coeffG <- nacopula:::coeffG

## ==== step (1): look at the a_k's, check if they can be evaluated ====


## Now, explore things seriously :
asN <- function(x, name=deparse(substitute(x))[1]) {
    names(x) <- paste(name, vapply(x, format, ""), sep="=")
    x
}
## use all methods for a set of alpha and d.vec
(meth <- eval(formals(coeffG)$method))
alpha <- c(.3, .5, .7, .8, .9, .99, .995)
d.vec <- c(5,seq(10, 90, by=5))

## *really* need the fixed sapply() {or R >= 2.13.x} !
## get the improved sapply():
if(getRversion() < "2.13")
    source(system.file("Rsource", "fixup-sapply.R", package="nacopula"))

ak.all <- sapply(asN(d.vec, "d"), function(d) {
    cat("\nd = ", d,"\n--------\n\n")
    sapply(asN(alpha), function(al) {
        cat("alpha = ", format(al), "\n")
        sapply(meth, coeffG, d=d, alpha=al, log=TRUE)
    }, simplify = "array")
}, simplify=FALSE)
## --> > 50 warnings

str(head(ak.all, 3))
ak.all$`d=20`[,,"alpha=0.99"]
chk1 <- function(ak.mat, tol = 1e-7) {
    stopifnot(is.matrix(ak.mat), (d <- nrow(ak.mat)) >= 2)
    n.meth <- ncol(ak.mat)
    med <- apply(ak.mat, 1, median, na.rm=TRUE)
    apply(ak.mat, 2, all.equal, target=med, tol=tol)
}

chk1(ak.all$`d=20`[,,"alpha=0.3"])
chk1(ak.all$`d=90`[,,"alpha=0.3"])

chk.all <-
    lapply(ak.all, function(ak.arr) apply(ak.arr, 3, chk1))

chk.all # quite interesting

## For d = 5,..85  this is fine (unless for large alpha (!) :

(a.k <- coeffG(100, 0.55, method = "horner"))
## => just works [but in the "extreme area", the numbers are not quite correct,
##    e.g., a.k[53] = 4.325e+83 and Maple says 4.627673570e83]

## conclusion: large alpha's [small theta's] cause the problems!!!
## ==========

## An example showing that for  "dsumSibuya" the problem is exactly *small* alphas:
plot (a.k.H <- coeffG(100, 0.01, method = "horner"), type = "l", lwd=3, log="y")
lines(a.k.J <- coeffG(100, 0.01, method = "dsumSibuya"), col=2, type ="o")
lines(a.k.s <- coeffG(100, 0.01, method = "sort"), col=3, type ="l")
lines(a.k.d <- coeffG(100, 0.01, method = "direct"),
      col=adjustcolor("blue"), type ="l", lwd=4)


set.seed(1)

n <- 50
d <- 100
tau <- 0.2
theta <- copGumbel@tauInv(tau)
alpha <- 1/theta

## animate this
library(animation)
library(lattice)

m <- 50 # frames
plot.list <- vector("list", m)
alpha.list <-  (1:m)/(m+1)
d <- 100
for(i in 1:m){
    coeffs <- coeffG(d, alpha.list[i], log=TRUE, method = "dsumSibuya")
    plot.list[[i]] <-
        xyplot(coeffs~1:d, type="l", xlim = c(-3,104), ylim = c(-303,374),
               xlab = "k", ylab = expression(log(a[k])), aspect = 1,
               main = substitute(expression(alpha == alpha.),
               list(alpha. = alpha.list[i])))
}
saveHTML(for(i in 1:m) print(plot.list[[i]]))

## conclusion: seems to be better for large alpha [seems to work for alpha >= 0.78,
##             including our test case by a whisker... not totally satisfactory so far]


plot(coeffs~1:100, type="l", xlim = c(-3,104), ylim = c(-303,374),
                             xlab = "k", ylab = expression(log(a[k])), aspect = 1,
                             main = substitute(expression(alpha == alpha.),
                             list(alpha. = alpha.list[i])))
