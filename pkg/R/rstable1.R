rstable1 <- function(n, alpha, beta, gamma = 1, delta = 0)
{
    ## Purpose: Random Numbers of Stable distribution in "1 - parametrization"
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, based on Diethelm Wuertz's code in fBasics

    stopifnot(0 < alpha, alpha <= 2,
              -1 <= beta, beta <= 1)

    ## if(pm == 1) ##  "S" aka "S1" or "1" Parametrization
    ## [[ if you drop this, you have the "S0" or "0" parametrization ]]
    delta <- delta + beta * gamma *
        ifelse(alpha == 1, (2/pi)*log(gamma), tan(pi*alpha/2))

    # Special case (a,b) = (1,0):
    if (all(alpha == 1) && all(beta == 0)) {
        res <- rcauchy(n)
    }
    else {
        ## Calculate uniform and exponential distributed random numbers:
        theta <- pi * (runif(n)-1/2)
        w <- rexp(n) ## -log(runif(n))

        ## FIXME: for alpha == 1, b.tan.pa <- Inf , but that needs different code
        ##        N.B.: Do this such that it vectorizes in alpha, too
        ## FIXME(2): ditto for    | alpha - 1 | << 1
        b.tan.pa <- beta*tan(pi*alpha/2)
        c. <- (1 + b.tan.pa^2)^(1/(2*alpha))
        th0 <- atan(b.tan.pa) / alpha
        a.tht <- alpha*(theta+th0)
        res <- c.* sin(a.tht) / cos(theta)^(1/alpha) *
            (cos(theta-a.tht)/w)^((1-alpha)/alpha)
        ## Using parametrization 0:
        res <- res - b.tan.pa
    }

    # Result:
    res * gamma + delta
}

