require(nacopula)

##' @title [m]inus Log Likelihood for Archimedean copula "fast version"
##' @param theta parameter (length 1 for our current families)
##' @param acop  Archimedean copula (of class "acopula")
##' @param u    data matrix n x d
##' @param n.MC NULL or positive integer for Monte-Carlo approximation
##' @return
##' @author Martin Maechler (Marius originally)
mLogL <- function(theta, acop, u, n.MC=NULL, ...) { # -(log-likelihood)
    -sum(acop@dacopula(u, theta, n.MC = n.MC, log = TRUE, ...))
}

##' @title
##' @param cop an 'outer_nacopula ' (currently with no children)
##' @param u n x d  data matrix
##' @param xlim x-range for curve() plotting
##' @param main title for curve()
##' @param XtrArgs a list of further arguments for mLogL()
##' @param ... further arguments for curve()
##' @return
##' @author Martin Maechler
curveLogL <- function(cop, u, xlim, main, XtrArgs=list(), ...) {
    unam <- deparse(substitute(u))
    stopifnot(is(cop, "outer_nacopula"), is.list(XtrArgs),
              (d <- ncol(u)) >= 2, d == dim(cop),
              length(cop@childCops) == 0# not yet *nested* A.copulas
              )
    acop <- cop@copula
    th. <- acop@theta
    if(missing(main)) {
        tau. <- cop@copula@tau(th.)
        main <- substitute("Neg. Log Lik."~ -italic(l)(theta, UU) ~ TXT ~~
			   FUN(theta['*'] == Th) %=>% tau['*'] == Tau,
			   list(UU = unam,
				TXT= sprintf("; n=%d, d=%d;  A.cop",
					     nrow(u), d),
				FUN = acop@name,
				Th = format(th.,digits=3),
				Tau = format(tau., digits=3)))
    }
    r <- curve(do.call(Vectorize(mLogL, "theta"), c(list(x, acop, u), XtrArgs)),
               xlim=xlim, main=main,
               xlab = expression(theta),
               ylab = substitute(- log(L(theta, u, ~~ COP)), list(COP=acop@name)),
               ...)
    axis(1, at = th., labels=expression(theta["*"]),
         lwd=2, col="dark gray", tck = -1/30)
    axis(1, at = paraOptInterval(u, acop@name),
         labels = FALSE, lwd = 2, col = 2, tck = 1/20)
    invisible(r)
}


n <- 200
d <- 100
tau <- 0.2
(theta <- copJoe@tauInv(tau))# 1.44381
(cop <- onacopulaL("Joe",list(theta,1:d)))

set.seed(1)
U1 <- rnacopula(n,cop)
enacopula(U1, cop, "mle") # 1.432885 --  fine
(th4 <- 1 + (1:4)/4)
mL.tr <- c(-3558.5, -3734.4, -3299.5, -2505.)
mLt1 <- sapply(th4, function(th) mLogL(th, cop@copula, U1))
mLt2 <- sapply(th4, function(th) mLogL(th, cop@copula, U1, method="maxSc"))
mLt3 <- sapply(th4, function(th) mLogL(th, cop@copula, U1, method="Jpoly"))
stopifnot(all.equal(mLt1, mL.tr, tol=5e-5),
          TRUE,## FIXME all.equal(mLt2, mL.tr, tol=5e-5),
          all.equal(mLt3, mL.tr, tol=5e-5))
##--> Funktion für Gesamtplot:
system.time(r1l  <- curveLogL(cop, U1, c(1, 2.5), X=list(method="logJpoly")))
system.time(r1J  <- curveLogL(cop, U1, c(1, 2.5), X=list(method="Jpoly"),
                              add=TRUE, col=adjustcolor("red", .4)))
system.time(r1m  <- curveLogL(cop, U1, c(1, 2.5), X=list(method="maxSc"),
                              add=TRUE, col=adjustcolor("blue",.5)))


U2 <- rnacopula(n,cop)
enacopula(U2, cop, "mle") # 1.4399  -- no warning any more
## the density for the *correct* parameter looks okay
summary(dnacopula(cop, U2))
## hmm:  max = 5.5e177
r2 <- curveLogL(cop, U2, c(1, 2.5))

mLogL(1.8, cop@copula, U2)# -4070.148 (was -Inf)

U3 <- rnacopula(n,cop)
enacopula(U3, cop, "mle") # 1.44957
r3 <- curveLogL(cop, U3, c(1, 2.5))

U4 <- rnacopula(n,cop)
enacopula(U4, cop, "mle") # 1.451929  was 2.351..  "completely wrong"
summary(dnacopula(cop, U4)) # ok (had one Inf)
r4 <- curveLogL(cop, U4, c(1, 2.5))

mLogL(2.2351, cop@copula, U4)
mLogL(1.5,    cop@copula, U4)
mLogL(1.2,    cop@copula, U4)

r4. <- curveLogL(cop, U4, c(1, 1.01))
r4. <- curveLogL(cop, U4, c(1, 1.0001))
r4. <- curveLogL(cop, U4, c(1, 1.000001))
##--> limit goes *VERY* steeply  up to .. probably 0

mLogL(1.2,    cop@copula, U4)
##--> theta 1.164 is about the boundary:
cop@copula@dacopula(U4[118,], theta=1.164, log = TRUE)
##  600.5926  (was Inf)
## now that we have  psiDabsJpoly(...., log=TRUE)

## Now a harder case:
n <- 200
d <- 150
tau <- 0.3
(theta <- copJoe@tauInv(tau))# 1.772
(cop <- onacopulaL("Joe",list(theta,1:d)))
U. <- rnacopula(n,cop)
enacopula(U., cop, "mle") # 1.776743
r. <- curveLogL(cop, U., c(1.1, 3))
## still looks very good
## and harder still

d <- 180
tau <- 0.4
(theta <- copJoe@tauInv(tau))# 2.219
(cop <- onacopulaL("Joe",list(theta,1:d)))
U. <- rnacopula(n,cop)
enacopula(U., cop, "mle") # 2.22666
r. <- curveLogL(cop, U., c(1.1, 4))
## still looks very good

1



###--------------------- The same for  Gumbel ---------------------------------
n <- 200
d <- 100
tau <- 0.2
(theta <- copGumbel@tauInv(tau))# 1.25
(cop <- onacopulaL("Gumbel",list(theta,1:d)))

set.seed(1)
U1 <- rnacopula(n,cop)
enacopula(U1, cop, "mle") # 1.53287 (not fine) -- 39 warnings
##--> Funktion für Gesamtplot:
r1 <- curveLogL(cop, U1, c(1, 2.1))

U2 <- rnacopula(n,cop)
r2 <- curveLogL(cop, U2, c(1, 2.1))

U3 <- rnacopula(n,cop)
r3 <- curveLogL(cop, U3, c(1, 2.1))

## all these similar: nowhere near minimum

mLogL(1.65, cop@copula, U1)# ok
mLogL(1.64, cop@copula, U1)# -> NaN
dd <- cop@copula@dacopula(U1, 1.64, log = TRUE)
dd ## viele  NaN ' s
cop@copula@dacopula(U1[10,], 1.64, log = TRUE)# einer jenen Fälle

debug(cop@copula@dacopula)
