require(nacopula)
options(warn = 1)

#### ==== Estimation and goodness-of-fit =========================================

source(system.file("Rsource", "estim-gof-fn.R", package="nacopula"))
## --> estimation.gof() etc

### ==== setup ====

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
    format(theta),", tau = ", format(tau),") ====\n\n",sep="")

if(getRversion() <= "2.13")
    source(system.file("Rsource", "fixup-sapply.R", package="nacopula"))

## For now, don't use the mde methods here, they take long
i.mde.meths <- grep("^mde", estMeth)
e.meths <- estMeth[-i.mde.meths]

RR <- sapply(e.meths, simplify="array", function(e)
         {
             n.MC <- if(e == "smle") 10000 else NULL
             sapply(gofMeth, simplify="array", function(g)
                    estimation.gof(n, d, simFamily, tau = tau, n.MC = n.MC,
                                   esti.method = e, gof.method = g))
         })

str(RR)
## Now print RR smartly (well, "to be improved"):
options(digits = 5)

## *Not* the times here:
RR[,c(1:2,4:5),,]

## but here:
apply(RR[,c("timeEstim","timeGOF"),,], c(4,1,2), mean)



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
if(!exists("doPlot")) doPlot <- TRUE

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

