#### Construct our "list" of supported Archimedean Copulas

### FIXME: Not "nice" that the names have to be used
###   probably: use *same* environment for  (tau, tauInv, paraConstr, nestConstr) ?

### ====Ali-Mikhail-Haq, see Nelsen (2007) p. 116, # 3====
copAMH <-
    new("ACopula", name = "AMH",
        ## generator:
        psi = function(t,theta) { (1-theta)/(exp(t)-theta) },
        ## generator inverse:
        psiInv = function(t,theta) { log((1-theta*(1-t))/t) },
        ## parameter constraint
 ## This is now *AUTOMAGICALLY* constructed (!) --> "initialize" in ./AllClass.R
 ## paraConstr = function(theta) {
 ##     length(theta) == 1 && 0 <= theta && theta < 1 },
        paraInterval = interval("[0,1)"),
        ## nesting constraint
        nestConstr = function(theta0,theta1) {
            copAMH@paraConstr(theta0) &&
            copAMH@paraConstr(theta1) && theta1 >= theta0
        },
        V0 = function(n,theta) { rgeom(n, 1-theta) + 1 },

        V01 = function(V0,theta0,theta1) { ## MM: FIXME?? do without sapply !
            variates <- sapply(V0, rgeom, prob = (1-theta1)/(1-theta0))
            sapply(variates,sum) + V0
        },

        tau = function(theta) {
            1 - 2*((1-theta)*(1-theta)*log(1-theta)+theta)/(3*theta*theta)
        },
        ## tau inverse
        tauInv = function(tau) {
            if(tau > 1/3) {
                stop("It is not possible for an Ali-Mikhail-Haq copula to attain such Kendall's tau")
            } else {
                r <- uniroot(function(t.) copAMH@tau(t.) - tau,
                             interval = c(1e-12, 1-1e-12))
                ## FIXME: check for convergence
                r$root
            }
        },

        ## lambda_l
        lTDC = function(theta) return(0*theta),
        ## lambda_l inverse
        lTDCInv = function(lambda) {
            if(any(lambda != 0))
                stop("Any parameter for an Ali-Mikhail-Haq copula gives a zero lower tail dependence coefficient")
            NA * lambda
        },

        ## lambda_u
        uTDC = function(theta) return(0*theta),
        uTDCInv = function(lambda) {
            if(any(lambda != 0))
                stop("Any parameter for an Ali-Mikhail-Haq copula gives a zero upper tail dependence coefficient")
            NA * lambda
        }
        )

##explicit check: stopifnot(validObject(copAMH))# ok

### ====Clayton, see Nelsen (2007) p. 116, # 1 but we use a slightly simpler form of the generator====

copClayton <-
    new("ACopula", name = "Clayton",
        ## generator and its inverse:
        psi = function(t,theta) { (1+t)^(-1/theta) },
        psiInv = function(t,theta) { t^(-theta) - 1 },
        paraInterval = interval("(0, Inf)"),
        ## nesting constraint
        nestConstr = function(theta0,theta1) {
            constrC(theta0) && constrC(theta1) && theta1 >= theta0
            ##----- FIXME
        },
        V0 = function(n,theta) { rgamma(n, shape = 1/theta) },
        V01 = function(V0,theta0,theta1) {
            alpha <- theta0/theta1
            retstable(alpha,V0,1)
        },
        tau = function(theta) { theta/(theta+2) },
        tauInv = function(tau) { 2*tau/(1-tau) },

        ## lambda_l and inverse:
        lTDC = function(theta) 2^(-1/theta),
        lTDCInv = function(lambda) -1/log2(lambda),

        ## lambda_u
        uTDC = function(theta) return(0*theta),
        uTDCInv = function(lambda) {
            if(any(lambda != 0))
                stop("Any parameter for a CLayton copula gives a zero upper tail dependence coefficient")
            NA * lambda
        })

stopifnot(validObject(copClayton))# ok

