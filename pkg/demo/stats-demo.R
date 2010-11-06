require(nacopula)

set.seed(1) # set seed

d <- 10
family.list <- c("AMH","Clayton","Frank","Gumbel","Joe")

## ==== Estimation =============================================================

## Demonstrates the fitting and goodness-of-fit capabilities for Archimedean
## copulas based on the provided families and dimension d.
estimation.gof <- function(n,family,tau,method = c("mle.tau.mean","mle.theta.mean",
                                    "mle.diag","smle","tau.tau.mean",
                                    "tau.theta.mean","dmle","beta")){

    ## generate data
    copFamily <- getAcop(family)
    theta <- copFamily@tauInv(tau)
    cop <- onacopulaL(family,list(theta,1:d))
    u <- rnacopula(n,cop)

    ## estimation and gof 
    num.families <- length(family.list)
    estimate <- numeric(num.families)
    gof <- numeric(num.families)
    for(k in 1:num.families){
        ## estimate the parameter with the provided method    
        cop.hat <- onacopulaL(family.list[k],list(,1:d))    
        estimate[k] <- enacopula(u,cop.hat,"mle.tau.mean",do.pseudo = FALSE)
	cop.hat@theta <- estimate[k]
	## apply a goodness-of-fit test to the estimated copula
        gof[k] <- gnacopula(u,cop.hat, do.pseudo = FALSE, ad.test = TRUE, 
                            method = "normal", MC = FALSE, bootstrap = FALSE)
    }
    
    ## result
    res <- cbind(estimate,gof)
    colnames(res) <- c("theta hat","p-value")
    row.names(res) <- family.list
    res

}

## Estimates for the different families based on n = 1000 d-dimensional vectors 
## of random variates from a Gumbel copula with Kendall's tau = 0.5
estimation.gof(1000,"Gumbel",0.5)
