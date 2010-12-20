require(nacopula)

##' @title [m]inus Log Likelihood for Archimedean copula "fast version"
##' @param theta parameter (length 1 for our current families)
##' @param acop  Archimedean copula (of class "acopula")
##' @param u    data matrix n x d
##' @param n.MC NULL or positive integer for Monte-Carlo approximation
##' @return
##' @author Martin Maechler (Marius originally)
mLogL <- function(theta, acop, u, n.MC=NULL) { # -(log-likelihood)
    -sum(acop@dacopula(u, theta, n.MC = n.MC, log = TRUE))
}

curveLogL <- function(cop, u, xlim, main, ...) {
    unam <- deparse(substitute(u))
    stopifnot(is(cop, "outer_nacopula"),
              (d <- ncol(u)) >= 2, d == dim(cop),
              length(cop@childCops) == 0# not yet *nested* A.copulas
              )
    acop <- cop@copula
    if(missing(main))
        main <- sprintf("Negative Log Likelihood  -l(., %s);  n=%d, d=%d;  A.copula '%s'",
                        unam, nrow(u), d, acop@name)
    r <- curve(Vectorize(mLogL, "theta")(x, acop, u), xlim=xlim,
               xlab = expression(theta),
               ylab = substitute(- log(L(theta, u, ~~ COP)), list(COP=acop@name)),
               main = main, ...)
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
##--> Funktion für Gesamtplot:
r1 <- curveLogL(cop, U1, c(1, 2.5))


U2 <- rnacopula(n,cop)
enacopula(U2, cop, "mle") # 1.4399  -- warning -- but that looks good to me (MM)
## the density for the *correct* parameter looks okay
summary(dnacopula(cop, U2))
## hmm:  max = 5.5e177
r2 <- curveLogL(cop, U2, c(1, 2.5))

mLogL(1.8, cop@copula, U2)# -Inf

U3 <- rnacopula(n,cop)
enacopula(U3, cop, "mle") # 1.44957
r3 <- curveLogL(cop, U3, c(1, 2.5))

U4 <- rnacopula(n,cop)
enacopula(U4, cop, "mle") # 2.351..  "completely wrong"
summary(dnacopula(cop, U4)) # one Inf
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
## Inf
##---> debug(..) clearly shows that the overflow to Inf happens
## inside   psiDabsJpoly(....)
##
## NB:  Instead of
        ## sum. <- psiDabsJpoly(arg, alpha, d)
        ## res[n01] <- (d - 1) * log(theta) + alpha * m.l.u + log(sum.) ....
## we need  psiDabsJpoly(...., log=TRUE)
##
##--> working via  sum-in-log-scale


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
