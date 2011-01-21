#### Is source()d from  demo(estimation.gof)  and also in some of the tests


## measuring user run time in milliseconds
utms <- function(x) 1000 * system.time(x)[[1]]
## formatting such times "uniformly":
f.tms <- function(x) paste(round(x),"ms") # with a space (sep = " ") !

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

    if(simFamily=="AMH"){
	n <- 200
	warning("FIXME in estimation.gof: changed n to 200 for AMH")
    }

    ## generate data
    copFamily <- getAcop(simFamily)
    theta <- copFamily@tauInv(tau)
    cop <- onacopulaL(simFamily, list(theta,1:d))
    if(verbose){
        r <- function(x) round(x,4) # for output
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
        ute[k] <- utms(est[k] <- enacopula(u, cop = cop.hat,
                                           method = esti.method,
                                           n.MC = n.MC, do.pseudo = FALSE))
        ## FIXME: test for "convergence" etc
        tau[k] <- cop.hat@copula@tau(est[k])
        if(verbose){
            cat("   theta hat      = ",r(est[k]),"\n",
                "   tau hat        = ",r(cop.hat@copula@tau(est[k])),"\n",
                ## The exact string 'Time ' must appear at the beginning of line
                ## for 'R CMD Rdiff' to *not* look at differences there:
                "Time estimation   = ",f.tms(ute[k]),"\n", sep="")
	}
        cop.hat@copula@theta <- est[k]
        ## apply a goodness-of-fit test to the estimated copula
        ## {{ use gnacopulatrafo() if you want the transformed u }}
        utg[k] <- utms(gof[k] <-
                       gnacopula(u, cop=cop.hat, estimation.method=esti.method,
                                 include.K=ncol(u)<=5, n.MC=n.MC, method=gof.method,
                                 do.pseudo = FALSE, do.pseudo.sim=FALSE, 
                                 verbose = FALSE))
        sig[k] <- if(gof[k] < 0.05) 1 else 0
        if(verbose){
            cat("   p-value        = ",r(gof[k]),"\n",
                "   < 0.05         = ",if(sig[k]) "TRUE" else "FALSE","\n",
                "Time G.O.Fit comp = ",f.tms(utg[k]),"\n\n", sep="")
	}
    }

    ## result
    names(est) <- checkFamilies
    cbind(theta_hat = est, tau_hat   = tau,
	  timeEstim = ute,
	  P_value   = gof, "< 0.05"  = sig,
	  timeGOF   = utg)
}

stopifnot(require(nacopula))
