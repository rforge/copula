require(nacopula)

## ==== Estimation =============================================================

## Demonstrates the fitting and goodness-of-fit capabilities for Archimedean
## copulas ...
## FIXME   : systematically try all 'method's -- not just one
## FIXME(2): systematically try all gnacopula-methods -- not just "normal"
estimation.gof <- function(n, d, dataFamily, tau,
                           checkFamilies = nacopula:::c_longNames,
                           MC=FALSE,
                           N = NULL,# < FIXME, something depending on family,theta,...
                           method = c("mle.tau.mean","mle.theta.mean",
                           "mle.diag","smle","tau.tau.mean",
                           "tau.theta.mean","dmle","beta"))
{
    ## generate data
    copFamily <- getAcop(dataFamily)
    theta <- copFamily@tauInv(tau)
    cop <- onacopulaL(dataFamily, list(theta,1:d))
    u <- rnacopula(n,cop)

    ## estimation and gof
    method <- match.arg(method)
    num.families <- length(checkFamilies)
    full.est <- as.list(estimate <- numeric(num.families))
    gof <- numeric(num.families)
    for(k in 1:num.families){
        ## estimate the parameter with the provided method
        cop.hat <- onacopulaL(checkFamilies[k],list(NA,1:d))
        est <- enacopula(u,cop.hat, method=method, do.pseudo = FALSE)
        ## FIXME: test for "convergence" etc
        full.est[[k]] <- est
        ## as long as enacopula() calls functions which call optimx() {{a mess!!}}
        ##             vvvvvv              vvvv {juck!}
        estimate[k] <- unlist(full.est[[k]]$par)
	cop.hat@copula@theta <- estimate[k]

	## apply a goodness-of-fit test to the estimated copula
        ## {{ use gnacopulatrafo() if you want the transformed u }}
        gof[k] <- gnacopula(u, cop.hat, do.pseudo = FALSE,
                            method = "normal", MC=MC, N=N, bootstrap = FALSE)
    }

    ## result
    res <- cbind(estimate,gof)
    colnames(res) <- c("theta hat","p-value")
    row.names(res) <- checkFamilies
    res
}

set.seed(1) # set seed

## Start with a "small" example:
estimation.gof(256, d = 3, "Gumbel", tau = 0.25, N = 2000)
## for now:
warnings()

## Estimates for the different families based on n = 1000 d-dimensional vectors
## of random variates from a Gumbel copula with Kendall's tau = 0.5
estimation.gof(1000, d=10, "Gumbel",
               tau = 0.5,
               ## as "AMH" cannot get a tau > 1/3 :
               checkFamilies = local({f <- nacopula:::c_longNames; f[f != "AMH"]}))

warnings()
