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

require(copula)

(doExtras <- copula:::doExtras())

if(!dev.interactive(orNone=TRUE)) pdf("asymCopula-ex.pdf")

### some constructions ###########################################################

## A Khoudraji-Clayton copula
kc <- khoudrajiCopula(copula2 = claytonCopula(6),
                      shapes = c(0.4, 0.95))
kc
contour(kc, dCopula, nlevels = 20, main = "dCopula(<khoudrajiBivCopula>)")

## check density: special case where we know the answer
kd2a <- khoudrajiCopula(copula1 = indepCopula(),
                        copula2 = claytonCopula(6),
                        shapes = c(1, 1))

kd2b <- khoudrajiCopula(copula1 = gumbelCopula(4),
                        copula2 = indepCopula(),
                        shapes = c(0, 0))

set.seed(12)
v <- matrix(runif(10), 5, 2)

stopifnot(all.equal(dCopula(v, kd2a), dCopula(v, claytonCopula(6))))
stopifnot(all.equal(dCopula(v, kd2b), dCopula(v, gumbelCopula(4))))


## True versus numerical derivatives
v <- matrix(runif(6), 3, 2)
max(abs(copula:::dCduCopulaNum(kc, v) - copula:::dCdu(kc, v)))
max(abs(copula:::dCdthetaCopulaNum(kc, v) - copula:::dCdtheta(kc, v)))

## tau, rho, lambda not supposed to work
assertError <- tools::assertError
assertError(tau(kc))
assertError(rho(kc))
assertError(iTau(kc, 0.5))
assertError(iRho(kc, 0.5))
assertError(lambda(kc))

## A Khoudraji-Clayton copula with one fixed shape parameter
kcf <- khoudrajiCopula(copula2 = claytonCopula(6),
                       shapes = fixParam(c(0.4, 0.95), c(FALSE, TRUE)))
kcf


## True versus numerical derivatives
v <- matrix(runif(6), 3, 2)
max(abs(copula:::dCduCopulaNum(kcf, v) - copula:::dCdu(kcf, v)))
max(abs(copula:::dCdthetaCopulaNum(kcf, v) - copula:::dCdtheta(kcf, v)))

## A Khoudraji-normal-Clayton copula
knc <- khoudrajiCopula(copula1 = normalCopula(-0.7),
                       copula2 = claytonCopula(6),
                       shapes = c(0.4, 0.95))
knc
contour(knc, dCopula, nlevels = 20, main = "dCopula(<khoudrajiBivCopula>)")

## True versus numerical derivatives
er.nc <- c(
dCdu = max(abs(copula:::dCduCopulaNum(knc, v) - copula:::dCdu(knc, v))),
dCdth= max(abs(copula:::dCdthetaCopulaNum(knc, v) - copula:::dCdtheta(knc, v))))
er.nc
##         dCdu        dCdth
## 0.0018449990 0.0001717451
stopifnot(abs(er.nc[["dCdu" ]]) < 0.004 ,
	  abs(er.nc[["dCdth"]]) < 0.0008)


## A Khoudraji-normal-Clayton copula with fixed params
kncf <- khoudrajiCopula(copula1 = normalCopula(fixParam(-0.7, TRUE)),
			copula2 = claytonCopula(6),
			shapes = fixParam(c(0.4, 0.95), c(FALSE, TRUE)))
kncf

## Test setTheta
kncf2 <- setTheta(kncf, value = c(0.5, 4, 0.2, 0.8))
kncf2

## True versus numerical derivatives
er.ncf <- c(
dCdu = max(abs(copula:::dCduCopulaNum(kncf, v) - copula:::dCdu(knc, v))),
dCdth= max(abs(copula:::dCdthetaCopulaNum(kncf, v) - copula:::dCdtheta(kncf, v))))
er.ncf
##         dCdu        dCdth
## 1.186238e-03 6.439294e-15
stopifnot(abs(er.ncf[["dCdu" ]]) < 0.004, # seen 2e-3
	  abs(er.ncf[["dCdth"]]) < 1e-12)

## A "nested" Khoudraji bivariate copula
kgkcf <- khoudrajiCopula(copula1 = gumbelCopula(3),
			 copula2 = kcf,
			 shapes = c(0.7, 0.25))
kgkcf
contour(kgkcf, dCopula, nlevels = 20, main = "dCopula(<khoudrajiBivCopula>)")
erN <- c(
dCdu = max(abs(copula:::dCduCopulaNum(kgkcf, v) - copula:::dCdu(kgkcf, v))),
dCdth= max(abs(copula:::dCdthetaCopulaNum(kgkcf, v) - copula:::dCdtheta(kgkcf, v))))
erN
##         dCdu        dCdth
## 1.340428e-03 8.121975e-14
stopifnot(abs(erN[["dCdu" ]]) < 3e-3,
	  abs(erN[["dCdth"]]) < 1e-12)

## A three dimensional Khoudraji-Clayton copula
kcd3 <- khoudrajiCopula(copula1 = indepCopula(dim=3),
                        copula2 = claytonCopula(6, dim=3),
                        shapes = c(0.4, 0.95, 0.95))
kcd3
class(kcd3) ## "khoudrajiExplicitCopula"
set.seed(1712)
n <- 300
u <- rCopula(n, kcd3)
v <- matrix(runif(15), 5, 3)
splom2(u)
(f.v <- dCopula(v, kcd3))
stopifnot(
    all.equal(f.v,
              c(0.23661418, 2.3571216, 0.13311822, 0.0051837845, 0.056044199),
              tol = 1e-7)# seen 1.69e-8
)
## check density: special case where we know the answer
kd3a <- khoudrajiCopula(copula1 = indepCopula(dim=3),
                        copula2 = claytonCopula(6, dim=3),
                        shapes = c(1, 1, 1))

kd3b <- khoudrajiCopula(copula1 = gumbelCopula(4, dim=3),
                        copula2 = indepCopula(dim=3),
                        shapes = c(0, 0, 0))
stopifnot(all.equal(dCopula(v, kd3a), dCopula(v, claytonCopula(6, dim=3))))
stopifnot(all.equal(dCopula(v, kd3b), dCopula(v, gumbelCopula(4, dim=3))))


## A four dimensional Khoudraji-Normal copula
knd4 <- khoudrajiCopula(copula1 = indepCopula(dim=4),
                        copula2 = normalCopula(.9, dim=4),
                        shapes = c(0.4, 0.95, 0.95, 0.95))
knd4
class(knd4) ## "khoudrajiCopula"
u <- rCopula(n, knd4)
splom2(u)
v <- matrix(runif(20), 5, 4)
assertError(dCopula(v, knd4)) ## should fail


## a nested explicit Khourdaji copula
kd3 <- khoudrajiCopula(copula1 = khoudrajiCopula(copula1 = gumbelCopula(4, dim=3),
                                                 copula2 = claytonCopula(6, dim=3),
                                                 shapes = c(.4, 0.95, 0.95)),
                       copula2 = frankCopula(10, dim=3),
                       shapes = c(.4, 0.95, 0.95))
kd3
stopifnot(identical(as.character(class(kd3)),
                    "khoudrajiExplicitCopula"))#  a khoudrajiExplicitCopula


## a nested Khoudraji copula whose construction should result
## in a khoudrajiCopula object, but not a khoudrajiExplicitCopula object

kd3kn <- khoudrajiCopula(copula1 = khoudrajiCopula(copula1 = gumbelCopula(4, dim=3),
                                                 copula2 = claytonCopula(6, dim=3),
                                                 shapes = c(.4, 0.95, 0.95)),
                         copula2 = normalCopula(.5, dim=3),
                         shapes = c(.4, 0.95, 0.95))
kd3kn
stopifnot(identical(as.character(class(kd3kn)),
                    "khoudrajiCopula"))  ## not a khoudrajiExplicitCopula


### fitting ###########################################################
n <- 300
set.seed(17)
u <- rCopula(n, kc)
plot(u)

if (doExtras)
{
    fk1 <- fitCopula(khoudrajiCopula(copula2 = claytonCopula()),
                     start = c(1.1, 0.5, 0.5), data = pobs(u),
                     optim.method = "Nelder-Mead", optim.control = list(trace = TRUE))
    fk1
    print( summary(fk1) )

    ## kcf : second shape parameter fixed to 0.95
    fkN <- fitCopula(kcf,
                     start = c(1.1, 0.5), data = pobs(u),
                     optim.method = "Nelder-Mead", optim.control = list(trace = TRUE))
    fkN
    print( summary(fkN) )
    print( cN <- confint(fkN) )
    fkB <- fitCopula(kcf,
                     start = c(1.1, 0.5), data = pobs(u),
                     optim.method = "BFGS", optim.control = list(trace = TRUE))
    print( summary(fkB) )
    stopifnot(
        all.equal(coef(fk1), c(c2.param = 5.42332, shape1 = 0.364467, shape2 = 0.868297),
                  tol = 1e-4), # seen 2.7e-7
        all.equal(coef(fkN), c(c2.param = 4.40525, shape1 = 0.389693), tol = 1e-4), # seen 3e-7
        all.equal(coef(fkN), coef(fkB), tol = 1e-3) # seen 1.19e-4
        ,
	all.equal(c(cN), c(2.246888, 0.2798682,
			   6.563614, 0.4995171), tol = 1e-5)# 5.4e-8
    )

    ## GOF example
    ## takes too long:
    ## gofCopula(kcf, x = u, start = c(1.1, 0.5), optim.method = "BFGS")
    set.seed(12)
    g.kcf <- gofCopula(kcf, x = u, start = c(1.1, 0.5), optim.method = "BFGS", sim = "mult")
    g.kcf
    stopifnot(inherits(g.kcf, "htest"),
              all.equal(g.kcf$p.value, 0.05544456, tol = 1e-3) # seen 8.2e-8
              )

    ## check size of mult GOF test briefly
    ## do1 <- function() {
    ##     u <- rCopula(n, kc)
    ##     gofCopula(kcf, x = u, start = c(1.1, 0.5), optim.method = "Nelder-Mead",
    ##               sim = "mult")$p.value
    ## }
    ## M <- 1000
    ## res <- replicate(M, do1())
    ## mean(res < 0.05)

    ## under the alternative
    set.seed(11); u <- rCopula(n, gumbelCopula(4))
    ## gofCopula(kcf, x = u, start = c(1.1, 0.5), optim.method = "BFGS")
    gA <- try(gofCopula(kcf, x = u, start = c(1.1, 0.5), optim.method = "BFGS", sim = "mult"))
    ## optim(*, "BFGS") sometimes not converging, but the above ex. does:
    stopifnot(inherits(gA, "htest"))
    ## and  "Nelder-Mead" does:
    set.seed(17)
    g2 <- gofCopula(kcf, x = u, start = c(1.1, 0.5), optim.method = "Nelder-M", sim = "mult")
    stopifnot(inherits(g2, "htest"),
	      all.equal(g2$p.value, 0.0004995005, tol = 1e-4))# seen 1e-9
}

### "Same" With rotated copula ###################################################
 ## GETR
if (doExtras)
{

    ## check density: special case where we know the answer
    kd2a <- khoudrajiCopula(copula1 = indepCopula(),
                            copula2 = rotCopula(claytonCopula(6)),
                            shapes = c(1, 1))
    kd2b <- khoudrajiCopula(copula1 = rotCopula(gumbelCopula(4)),
                            copula2 = indepCopula(),
                            shapes = c(0, 0))
    v <- matrix(runif(10), 5, 2)
    stopifnot(all.equal(dCopula(v, kd2a), dCopula(v, rotCopula(claytonCopula(6)))))
    stopifnot(all.equal(dCopula(v, kd2b), dCopula(v, rotCopula(gumbelCopula(4)))))

    ## same with a different flip
    T.F. <- c(TRUE, FALSE)
    kd2a <- khoudrajiCopula(copula1 = indepCopula(),
                            copula2 = rotCopula(claytonCopula(6), flip = T.F.),
                            shapes = c(1, 1))
    kd2b <- khoudrajiCopula(copula1 = rotCopula(gumbelCopula(4), flip = T.F.),
                            copula2 = indepCopula(),
                            shapes = c(0, 0))
    v <- matrix(runif(10), 5, 2)
    stopifnot(all.equal(dCopula(v, kd2a),
                        dCopula(v, rotCopula(claytonCopula(6), flip = T.F.))))
    stopifnot(all.equal(dCopula(v, kd2b),
                        dCopula(v, rotCopula(gumbelCopula(4), flip = T.F.))))

    ## other basic checks
    F.T. <- c(FALSE, TRUE)
    krc <- khoudrajiCopula(copula1 = rotCopula(claytonCopula(6), flip = F.T.),
                           shapes = c(0, 0))
    contour(krc, dCopula, nlevels = 20, main = "dCopula(<khoudrajiBivCopula>)")
    rkc <- rotCopula(khoudrajiCopula(copula2 = claytonCopula(6),
                     shapes = c(1, 1)), flip = F.T.)
    contour(rkc, dCopula, nlevels = 20, main = "dCopula(<khoudrajiBivCopula>)")
    stopifnot(all.equal(dCopula(v, krc), dCopula(v, rkc)),
	      all.equal(pCopula(v, krc), pCopula(v, rkc)))


    ## Fitting
    ## A Khoudraji copula constructed from the survival Clayton
    krc <- khoudrajiCopula(copula2 = rotCopula(claytonCopula(6)),
                           shapes = c(0.4, 0.95))
    krc
    contour(krc, dCopula, nlevels = 20, main = "dCopula(<khoudrajiBivCopula>)")
    n <- 300
    set.seed(17)
    u <- rCopula(n, krc)
    plot(u)

    fk1 <- fitCopula(khoudrajiCopula(copula2 = rotCopula(claytonCopula())),
                     start = c(1.1, 0.5, 0.5), data = pobs(u),
                     optim.method = "Nelder-Mead", optim.control = list(trace = TRUE))
    print( summary(fk1) )


    ## second shape parameter fixed to 0.95
    fixedParam(krc) <- c(FALSE, FALSE, TRUE)
    fkN <- fitCopula(krc,
                     start = c(1.1, 0.5), data = pobs(u),
                     optim.method = "Nelder-Mead", optim.control = list(trace = TRUE))
    print( summary(fkN) )
    fkB <- fitCopula(krc,
                     start = c(1.1, 0.5), data = pobs(u),
                     optim.method = "BFGS", optim.control = list(trace = TRUE))
    print( summary(fkB) )
    stopifnot(
        all.equal(coef(fk1), c(c2.param = 4.328559, shape1 = 0.5109215, shape2 = 0.9653532),
                  tol = 1e-4) # seen 9.2e-8
       , all.equal(coef(fkN), c(c2.param = 4.372715, shape1 = 0.5042553), tol = 1e-4)# seen 3e-8
       , all.equal(coef(fkN), coef(fkB), tol = 4e-4) # seen 4.11e-5
    )

    ## GOF example
    ## takes too long:
    ## gofCopula(krc, x = u, start = c(1.1, 0.5), optim.method = "BFGS")
    set.seed(12)
    g.krc <- gofCopula(krc, x = u, start = c(1.1, 0.5), optim.method = "BFGS", sim = "mult")
    g.krc
    stopifnot(inherits(g.krc, "htest"),
              all.equal(g.krc$p.value, 0.32317682, tol = 1e-4) # seen 8.8e-9
              )

    ## check size of mult GOF test briefly
    ## do1 <- function() {
    ##     u <- rCopula(n, kc)
    ##     gofCopula(krc, x = u, start = c(1.1, 0.5), optim.method = "Nelder-Mead",
    ##               sim = "mult")$p.value
    ## }
    ## M <- 1000
    ## res <- replicate(M, do1())
    ## mean(res < 0.05)

    ## under the alternative
    set.seed(11); u <- rCopula(n, gumbelCopula(4))
    ## gofCopula(krc, x = u, start = c(1.1, 0.5), optim.method = "BFGS")
    gA <- try(gofCopula(krc, x = u, start = c(1.1, 0.5), optim.method = "BFGS", sim = "mult"))
    ## does not converge
    set.seed(17)
    g2 <- gofCopula(krc, x = u, start = c(1.1, 0.5), optim.method = "Nelder-M", sim = "mult")
    stopifnot(inherits(g2, "htest"),
              all.equal(g2$p.value, 0.0014985015, tol = 1e-4)# seen 1e-9
              )
}

if (FALSE) ## too time consuming .. NB 'estimate.variance'
{
    ##############################################################################
    ## Testing density / fitting

    fitest <- function(n, cop, N)
    {
        stopifnot(N >= 2)
        st1 <- function(fc) cbind(estim = coef(fc), sd.err = sqrt(diag(fc@var.est)))
        do1 <- function(stat.only = TRUE) {
            u <- rCopula(n, cop)
            fit <- fitCopula(cop, start = c(1.1, rep(0.5, dim(cop))), data = pobs(u),
                             optim.method = "Nelder-Mead", estimate.variance = TRUE)
            if(stat.only) st1(fit) else fit
        }
        ## Fit one:
        re1 <- do1(stat.only=FALSE)
        print(summary(re1))

        ## bias and stderr
        rr <- replicate(N-1, do1())
        ## here and below, I'd like abind():
        res <- array(c(st1(re1), rr), dim = dim(rr) + c(0L,0:1), dimnames=dimnames(rr))
        m <- apply(res, 1:2, mean)
        s <- apply(res, 1:2, sd)
        array(c(t(m),t(s)), dim = c(2, dim(m)), dimnames = c(list(c("Mean","SD")), dimnames(m)))
    }


    n <- 128 # want speed
    ##############################################################################
    ## Testing bivariate rotated / khoudraji

    cop <- rotCopula(khoudrajiCopula(copula1 = indepCopula(),
                                     copula2 = claytonCopula(6),
                                     shapes = c(0.4, 0.95)))
    set.seed(1251)
    system.time(fi2KC <- fitest(n=n, cop=cop, N=4))
    print(      fi2KC)

    ##############################################################################
    ## Testing bivariate khoudraji / rotated

    set.seed(1251)
    cop <- khoudrajiCopula(copula1 = indepCopula(),
                           copula2 = rotCopula(gumbelCopula(6)),
                           shapes = c(0.4, 0.95))
    system.time(fi2KrC <- fitest(n=n, cop=cop, N=4))
    print(      fi2KrC)


    ##############################################################################
    ## Testing higher-dimensional density / fitting

    ## A three dimensional Khoudraji-Clayton copula
    cop <- khoudrajiCopula(copula1 = indepCopula(dim=3),
                           copula2 = claytonCopula(6, dim=3),
                           shapes = c(0.4, 0.95, 0.95))
    system.time(fi3KC <- fitest(n=n, cop=cop, N=4))
    print(fi3KC)

    ## A three dimensional Khoudraji-Gumbel copula
    cop <- khoudrajiCopula(copula1 = indepCopula(dim=3),
                           copula2 = gumbelCopula(6, dim=3),
                           shapes = c(0.4, 0.95, 0.95))
    system.time(fi3KG <- fitest(n=n, cop=cop, N=4))
    print(fi3KG)
}

## All 'copula' subclasses
