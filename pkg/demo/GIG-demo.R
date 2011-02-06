library(nacopula) 
options(warn = 1) 

#### ==== Demo for working with the two-parameter Archimedean GIG copula =======

source(system.file("Rsource", "GIG.R", package="nacopula"))

## ==== setup ==================================================================

## set seed
set.seed(1) 

## ==== specify parameters ====

d <- 5 # dimension
tau <- 0.5 # specify Kendall's tau
## specify the theta vector such that Kendall's tau equals the given tau
theta <- tauInv.GIG(tau, theta=c(0.01, NA), interval=c(1e-12,100))

## initial interval bounds for optimization
h <- c(0.05, 0.05) # for initial interval for theta_1 and theta_2
I <- matrix(c(theta[1]-h[1], theta[1]+h[1], theta[2]-h[2], theta[2]+h[2]), ncol=2) # initial intervals for theta_1 and theta_2

## ==== sample the copula ====

n <- 1000
U <- rnacopula.GIG(n, d, theta)
tau.hat <- cor(U, method="kendall")
tau.hats <- tau.hat[upper.tri(tau.hat)]
stopifnot(all.equal(rep(tau, d*(d-1)/2), tau.hats, tol=0.03)) # check tau

## plot the sample
pairs(U, gap=0)

## ==== estimate the parameters ====

n <- 100 # choose small sample size due to run time [depending on psi^{-1}]
U <- rnacopula.GIG(n, d, theta)

## choose random initial point in a reasonable interval
start <- c(runif(1, min=I[1,1], max=I[1,2]), runif(1, min=I[2,1], max=I[2,2]))

## call optimizer for MLE
system.time(res.MLE <- optimx(par=start, 
                              fn=function(x) -sum(dacopula.GIG(U, x, n.MC=0, log=TRUE)), 
                              lower=c(I[1,1],I[2,1]), upper=c(I[1,2],I[2,2]), 
                              method="bobyqa"))
res.MLE

## call optimizer for SMLE
system.time(res.SMLE <- optimx(par=start, 
                               fn=function(x) -sum(dacopula.GIG(U, x, n.MC=10000, log=TRUE)), 
                               lower=c(I[1,1],I[2,1]), upper=c(I[1,2],I[2,2]), 
                               method="bobyqa"))
res.SMLE

## note: run time is mainly determined by the evaluations of psiInv.GIG

## ==== plot Kendall's tau as a function in theta for different nu's ====

th1 <- c(0, 0.001, 0.5, 1, 5, 10)
cols <- colorRampPalette(c("red", "orange", "darkgreen", "turquoise", "blue"), 
                         space="Lab")(length(th1))
for(i in seq_along(th1))
    curve(tau.GIG(cbind(th1[i],x)), 1e-12, 2, 
          main="Kendall's tau for the GIG family", ylim=c(0,1),
          xlab=expression(theta), ylab=expression(tau(nu,theta)), add=(i>1), 
          lwd=1.4, col=cols[i])
label <- as.expression(lapply(1:length(th1), function(i) substitute(nu==nu., list(nu.=th1[i]))))
legend("topright", label, bty="n", lwd=1.4, col=cols)
## conclusion: - largest range of tau's for theta_1 = 0
##             - not possible to evaluate for theta_1 < 0
##             - tau.GIG is numerically critical for theta_1 > 0 close to 0             
