## ##' @title Bivariate Conditional Copulas C(u2|u1) of u2 given u1
## ##' @param u2 data vector (numeric(n))
## ##' @param u1 data vector (numeric(n))
## ##' @param family copula family
## ##' @param theta parameter (for ACs; for elliptical copulas, its rho)
## ##' @param ... additional args (e.g., df for t copulas)
## ##' @return C(u2|u1)
## ##' @author Marius Hofert and Martin Maechler
## ##' Note: used in Hofert and Maechler (2013)
## ccop2 <- function(u2, u1, family, theta, ...) {
##     stopifnot(length(u1)==length(u2))
##     switch(copFamilyClass(family),
##            "ellipCopula"={
##                switch(family,
##                       "normal"={
##                           pnorm((qnorm(u2)-theta*qnorm(u1))/sqrt(1-theta^2))
##                       },
##                       "t"={
##                           stopifnot(hasArg(df))
##                           df <- list(...)$df
##                           qt.u1 <- qt(u1, df=df)
##                           mu <- theta * qt.u1
##                           sigma <- sqrt((df+qt.u1^2)*(1-theta^2)/(df+1))
##                           pt((qt(u2, df=df)-mu)/sigma, df=df+1)
##                       },
##                       stop("not yet supported family"))
##            },
##            "outer_nacopula"={
##                ## cacopula(u, cop=onacopulaL(family, list(theta, 1:2)))
##                cop <- getAcop(family)
##                u <- cbind(u1, u2)
##                psiI <- cop@iPsi(u, theta=theta)
## 	       exp(cop@absdPsi(rowSums(psiI), theta=theta, log=TRUE) -
## 		   cop@absdPsi(psiI[,1], theta=theta, log=TRUE))
##            },
##            stop("family ", family, " not yet supported"))
## }
