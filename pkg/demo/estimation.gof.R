require(nacopula)
options(warn = 1)

## ==== Estimation and goodness-of-fit =========================================

r <- function(x) round(x,4) # for output

##' Demonstrates the fitting and goodness-of-fit capabilities for Archimedean
##' copulas
##'
##' @title Fitting and Goodness-Of-Fit for Archimedean copulas
##' @param n sample size
##' @param d dimension
##' @param simFamily Archimedean family to be sampled
##' @param tau degree of dependence of the sampled family in terms of Kendall's tau
##' @param n.MC if provided (and not NULL) and if it make sense for the
##' chosen method, Monte Carlo is used with sample size equal to n.MC
##' @param esti.method estimation method (see enacopula)
##' @param gof.method  goodness-of-fit transformation (see gnacopula)
##' @param checkFamilies vector of Archimedean families to be used for gof
##' @param verbose
##' @return a numeric matrix ...
##' @author Marius Hofert and Martin Maechler
estimation.gof <- function(n, d, simFamily, tau, n.MC,
                           esti.method = eval(formals(enacopula)$method),
                           gof.method = eval(formals(gnacopula)$method),
                           checkFamilies = nacopula:::c_longNames, verbose = TRUE)
{

    ## generate data
    copFamily <- getAcop(simFamily)
    theta <- copFamily@tauInv(tau)
    cop <- onacopulaL(simFamily, list(theta,1:d))
    if(verbose){
	cat("\n\n## ==== Output for esti.method = \"",esti.method,
            "\" and gof.method = \"",gof.method,"\" ====\n\n",sep="")
    }
    u <- rnacopula(n,cop)

    ## estimation and gof
    esti.method <- match.arg(esti.method)
    gof.method <- match.arg(gof.method)
    n.cf <- length(checkFamilies)
    tau <- gof <- sig <- ute <- utg <- est <- numeric(n.cf)
    for(k in 1:n.cf) {
        ## estimate the parameter with the provided method
        cop.hat <- onacopulaL(checkFamilies[k],list(NA,1:d))
        if(verbose) cat("Estimation and GOF for ",checkFamilies[k],":\n\n",sep="")
        ute[k] <- system.time(est[k] <- enacopula(u, cop = cop.hat,
                                                   method = esti.method,
                                                   n.MC = n.MC, do.pseudo = FALSE))[1]
        ## FIXME: test for "convergence" etc
        tau[k] <- cop.hat@copula@tau(est[k])
        if(verbose){
            cat("   theta hat      = ",r(est[k]),"\n",
                "   tau hat        = ",r(cop.hat@copula@tau(est[k])),"\n",
            ## The exact string 'Time ' must appear at the beginning of line
            ## for 'R CMD Rdiff' to *not* look at differences there:
                "Time estimation   = ",round(1000*ute[k]),"ms\n",sep="")
	}
        cop.hat@copula@theta <- est[k]
        ## apply a goodness-of-fit test to the estimated copula
        ## {{ use gnacopulatrafo() if you want the transformed u }}
        utg[k] <-
            system.time(gof[k] <-
                        gnacopula(u, cop = cop.hat, bootstrap = FALSE,
                                  method = gof.method,
                                  estimation.method = esti.method,
                                  n.MC = n.MC, do.pseudo = FALSE,
                                  verbose = FALSE))[1]
        sig[k] <- if(gof[k] < 0.05) 1 else 0
        if(verbose){
            cat("   p-value        = ",r(gof[k]),"\n",
                "   < 0.05         = ",if(sig[k]) "TRUE" else "FALSE","\n",
                "Time G.O.Fit comp = ",round(1000*utg[k]),"ms\n\n",sep="")
	}
    }

    ## result
    names(est) <- checkFamilies
    cbind(theta_hat = est, tau_hat   = tau,
	  timeEstim = ute,
	  P_value   = gof, "< 0.05"  = sig,
	  timeGOF   = utg)
}


## ==== setup ====

## Use all available estimation and GOF methods:
(estMeth <- eval(formals(enacopula)$method))
(gofMeth <- eval(formals(gnacopula)$method))

set.seed(1) # set seed

n <- 256 # sample size
d <- 5 # dimension
tau <- 0.2 # Kendall's tau

## ==== apply all procedures (to data from AMH) ================================

simFamily <- "AMH"
cop <- getAcop(simFamily)
theta <- cop@tauInv(tau) # true parameter

## start the loop
cat("\n## ==== data from ",simFamily," (n = ",n,", d = ",d,", theta = ",
    r(theta),", tau = ",r(tau),") ====\n\n",sep="")
## sapply(*, simplify=FALSE):
## sapply() cannot (yet!) produce ## higher-level arrays
## Martin's patched version
sapply <- function(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE, ARRAY = FALSE)
{
    FUN <- match.fun(FUN)
    answer <- lapply(X, FUN, ...)
    if(USE.NAMES && is.character(X) && is.null(names(answer)))
	names(answer) <- X
    if(simplify && length(answer) &&
       length(common.len <- unique(unlist(lapply(answer, length)))) == 1L) {
	if(common.len == 1L)
	    unlist(answer, recursive = FALSE)
	else if(common.len > 1L) {
	    r <- as.vector(unlist(answer, recursive = FALSE))
	    if(ARRAY && length(c.dim <- unique(lapply(answer, dim))) == 1 &&
	       is.numeric(c.dim <- c.dim[[1L]]) &&
	       prod(d <- c(c.dim, length(X))) == length(r)) {

		iN1 <- is.null(n1 <- dimnames(answer[[1L]]))
                n2 <- names(answer)
		dnam <-
		    if(!(iN1 && is.null(n2)))
			c(if(iN1) rep.int(list(n1), length(c.dim)) else n1,
                          list(n2)) ## else NULL
		array(r, dim = d, dimnames = dnam)

	    } else if(prod(d <- c(common.len, length(X))) == length(r))

		array(r, dim = d,
		      dimnames= if(!(is.null(n1 <- names(answer[[1L]])) &
				     is.null(n2 <- names(answer)))) list(n1,n2))
	    else answer
	}
	else answer
    } else answer
}


RR <- sapply(estMeth, ARRAY=TRUE, function(e)
         {
             n.MC <- if(e == "smle") 10000 else NULL
             sapply(gofMeth, ARRAY=TRUE, function(g)
                    estimation.gof(n, d, simFamily, tau = tau, n.MC = n.MC,
                                   esti.method = e, gof.method = g))
         })

str(RR)
## Now print RR smartly (well, "to be improved"):
options(digits = 5)

## *Not* the times here:
RR[,c(1:2,4:5),,]

## but here:
apply(1000 * RR[,c("timeEstim","timeGOF"),,],
      c(4,1,2), mean)



## ==== MLE estimation by hand (for debugging purposes) ====

## generate the data
simFamily <- getAcop("AMH") # choose your desired family
cop <- onacopulaL(simFamily, list(theta <- simFamily@tauInv(tau),1:d))
u <- rnacopula(n,cop)

## estimate the copula
copFamily <- "Joe" # family to be estimated
cop.hat <- onacopulaL(copFamily,list(NA,1:d))
mLogL <- function(theta,cop.hat,u){
    cop.hat@copula@theta <- theta
    -sum(dnacopula(cop.hat, u, log = TRUE))
}
(est <- optimize(mLogL, interval = paraOptInterval(u,copFamily), cop.hat = cop.hat,
                 u = u))

## evaluate the density at u for a specified parameter
theta <- 14
cop.hat@copula@theta <- theta
(log.dens <- dnacopula(cop.hat, u, log = TRUE))
-sum(log.dens)

## ==== Plots ==================================================================

if(!dev.interactive())## e.g. when run as part of R CMD check
    pdf("demo_est-gof.pdf")
doPlot <- TRUE

## ==== setup for plots ====

t01 <- (0:256)/256
copA <- setTheta(copAMH, 0.7135001)
copC <- setTheta(copClayton, 0.5)
copF <- setTheta(copFrank, 1.860884)
copG <- setTheta(copGumbel, 1.25)
copJ <- setTheta(copJoe, 1.25)

cols <- c("black","orange3","red3","darkgreen","blue") # no very light ones
## TODO: work with *list* of copulas, and extra names *and* theta's (!)
##       to be also put in labels
labs <- c("AMH","Clayton","Frank","Gumbel","Joe")

## ==== plots of the densities of the diagonals ====

d <- 5

dDmat <- cbind(dDiag.A = dDiag(t01,copA,d),
               dDiag.C = dDiag(t01,copC,d),
               dDiag.F = dDiag(t01,copF,d),
               dDiag.G = dDiag(t01,copG,d),
               dDiag.J = dDiag(t01,copJ,d))

if(doPlot) {
    matplot(t01, dDmat, type="l", col=cols, xlab="t", ylab="dDiag(t)")
    legend("bottomright", legend=labs, lty=1:5, col=cols, bty="n")
    ## and in log-log scale:
    matplot(t01, dDmat, type="l", col=cols, xlab="t",
            log = "xy", main = "dDiag(t) [log-log scale]")
    legend("bottomright", legend=labs, lty=1:5, col=cols, bty="n")
}

## ==== plots of the Kendall distribution functions ====

d <- 10
Kmat <- cbind(K.A = K(t01,copA,d),
              K.C = K(t01,copC,d),
              K.F = K(t01,copF,d),
              K.G = K(t01,copG,d),
              K.J = K(t01,copJ,d))
head(mm <- cbind(t = t01, Kmat))
tail(mm)

dK <- apply(Kmat, 2, diff)
summary(dK)
## NOTE:  AMH and Clayton have some (very slightly) negative values
## ----    <==>  K() is not increasing there (near 1)
## MM: this is  "unavoidable" because of the numerics behind ...

if(doPlot) {
    matplot(t01, Kmat, type="l", col=cols, xlab="t", ylab="K(t)")
    legend("bottomright", legend=labs, lty=1:5, col=cols, bty="n")
    ## and in log-log scale:
    matplot(t01, Kmat, type="l", col=cols, xlab="t",
            log = "xy", main = "K(t) [log-log scale]")
    legend("bottomright", legend=labs, lty=1:5, col=cols, bty="n")
}

