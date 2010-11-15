require(nacopula)
options(warn = 1)

## ==== Estimation and goodness-of-fit =========================================

r <- function(x) round(x,4) # for output

## Demonstrates the fitting and goodness-of-fit capabilities for Archimedean
## copulas
## n sample size
## d dimension
## simFamily Archimedean family to be sampled
## tau degree of dependence of the sampled family in terms of Kendall's tau
## MC if provided (and not NULL) and if it make sense for the chosen method,
##    Monte Carlo is used with sample size equal to MC
## estimation.method estimation method (s. enacopula)
## gof.method goodness-of-fit transformation (s. gnacopula)
## checkFamilies vector of Archimedean families to be used for gof
estimation.gof <- function(n, d, simFamily, tau,
                           estimation.method = c(
                           "mle",
                           "smle",
                           "tau.tau.mean",
                           "tau.theta.mean",
                           "dmle",
                           "beta"
                           ),
                           MC,
			   gof.method = c(
			   "log",
                           "normal"
                           ),
                           checkFamilies = nacopula:::c_longNames, verbose = TRUE)
{

    ## generate data
    copFamily <- getAcop(simFamily)
    theta <- copFamily@tauInv(tau)
    cop <- onacopulaL(simFamily, list(theta,1:d))
    if(verbose){
	cat("\n\n## ==== Output for estimation.method = \"",estimation.method,
            "\" and gof.method = \"",gof.method,"\" ====\n\n",sep="")
    }
    u <- rnacopula(n,cop)

    ## estimation and gof
    estimation.method <- match.arg(estimation.method)
    gof.method <- match.arg(gof.method)
    num.checkFamilies <- length(checkFamilies)
    est <- numeric(num.checkFamilies)
    tau <- est
    gof <- est
    sig <- est
    ut <- matrix(,nrow=num.checkFamilies,ncol=2)
    for(k in 1:num.checkFamilies){
        ## estimate the parameter with the provided method
        cop.hat <- onacopulaL(checkFamilies[k],list(NA,1:d))
        if(verbose) cat("Estimation and GOF for ",checkFamilies[k],":\n\n",sep="")
        ut[k,1] <- system.time(est[k] <- enacopula(u, cop = cop.hat,
                                                   method = estimation.method,
                                                   MC = MC, do.pseudo = FALSE))[1]
        ## FIXME: test for "convergence" etc
        tau[k] <- cop.hat@copula@tau(est[k])
        if(verbose){
            cat("   theta hat      = ",r(est[k]),"\n",sep="")
            cat("   tau hat        = ",r(cop.hat@copula@tau(est[k])),"\n",sep="")
            cat("   user time est. = ",r(ut[k,1]),"s\n",sep="")
	}
        cop.hat@copula@theta <- est[k]
        ## apply a goodness-of-fit test to the estimated copula
        ## {{ use gnacopulatrafo() if you want the transformed u }}
        ut[k,2] <- system.time(gof[k] <- gnacopula(u, cop = cop.hat, bootstrap = FALSE,
                                                   method = gof.method,
                                                   estimation.method = estimation.method,
                                                   MC = MC, do.pseudo = FALSE,
                                                   verbose = FALSE))[1]
        sig[k] <- if(gof[k] < 0.05) 1 else 0
        if(verbose){
            cat("   p-value        = ",r(gof[k]),"\n",sep="")
            cat("   < 0.05         = ",if(sig[k]) "TRUE" else "FALSE","\n",sep="")
            cat("   user time gof  = ",r(ut[k,2]),"s\n\n",sep="")
	}
    }

    ## result
    res <- cbind(est,tau,ut[,1],gof,sig,ut[,2])
    colnames(res) <- c("theta hat","tau hat","user time est.","p-value","< 0.05",
                       "user time gof")
    row.names(res) <- checkFamilies
    res
}

## ==== setup ====

set.seed(1) # set seed

n <- 512 # sample size
d <- 5 # dimension
tau <- 0.25 # Kendall's tau
estMeth <- c("mle","smle","tau.tau.mean","tau.theta.mean","dmle","beta")
gofMeth <- c("log","normal")

## ==== apply all procedures (to data from AMH) ================================

simFamily <- "AMH"
cop <- getAcop(simFamily)
theta <- cop@tauInv(tau) # true parameter

## start the loop
cat("\n## ==== data from ",simFamily," (n = ",n,", d = ",d,", theta = ",
    r(theta),", tau = ",r(tau),") ====\n\n",sep="")
for(e in estMeth){
    MC <- if(e == "smle") 10000 else NULL
    for(g in gofMeth){
        res <- estimation.gof(n, d, simFamily, tau, estimation.method = e,
                              MC, gof.method = g)
        cat("Results:\n\n")
        print(r(res))
    }
}

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

