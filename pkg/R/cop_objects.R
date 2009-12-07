#### Construct our "list" of supported Archimedean Copulas

### FIXME: Not "nice" that the names have to be used
###   probably: use *same* environment for  (tau, tauInv, paraConstr, nestConstr) ?

### ====Ali-Mikhail-Haq, see Nelsen (2007) p. 116, # 3====
copAMH <-
    new("ACopula", name = "AMH",
        ## generator
        psi = function(t,theta) { (1-theta)/(exp(t)-theta) },
        psiInv = function(t,theta) { log((1-theta*(1-t))/t) },
        ## parameter interval
        paraInterval = interval("[0,1)"),
        ## nesting constraint
        nestConstr = function(theta0,theta1) {
            copAMH@paraConstr(theta0) &&
            copAMH@paraConstr(theta1) && theta1 >= theta0
        },
        ## V0 and V01
        V0 = function(n,theta) { rgeom(n, 1-theta) + 1 },
        V01 = function(V0,theta0,theta1) { ## MM: FIXME?? do without sapply !
            ##MH: this should be a V0-fold sum of rgeom(prob) random variates
            ##MH: for 1 V0, this should be sum(rgeom(V0,prob=(1-theta1)/(1-theta0)))
            variates <- sapply(V0, rgeom, prob = (1-theta1)/(1-theta0))
            sapply(variates,sum) + V0
        },
        ## Kendall's tau
        tau = function(theta) {
            1 - 2*((1-theta)*(1-theta)*log(1-theta)+theta)/(3*theta*theta)
        },
        tauInv = function(tau) {
                if(any(tau>1/3)){
                  stop("It is not possible for an AMH copula to attain a Kendall's tau larger than 1/3")
                }
                sapply(tau,function(tau){
                                    r <- uniroot(function(th) copAMH@tau(th) - tau,
                                                 interval = c(1e-12, 1-1e-12))
                                    r$root## FIXME: check for convergence
                                    })
        },     
        ## lower tail dependence coefficient lambda_l
        lTDC = function(theta) { 0*theta },
        lTDCInv = function(lambda) {
            if(any(lambda != 0))
                stop("Any parameter for an Ali-Mikhail-Haq copula gives zero lower tail dependence coefficient")
            NA * lambda
        },
        ## upper tail dependence coefficient lambda_u
        uTDC = function(theta) { 0*theta },
        uTDCInv = function(lambda) {
            if(any(lambda != 0))
                stop("Any parameter for an Ali-Mikhail-Haq copula gives zero upper tail dependence coefficient")
            NA * lambda
        })

stopifnot(validObject(copAMH))# ok

### === Clayton, see Nelsen (2007) p. 116, #1
### 	but we use a slightly simpler form of the generator====
copClayton <-
    new("ACopula", name = "Clayton",
        ## generator
        psi = function(t,theta) { (1+t)^(-1/theta) },
        psiInv = function(t,theta) { t^(-theta) - 1 },
        ## parameter interval
        paraInterval = interval("(0,Inf)"),
        ## nesting constraint
        nestConstr = function(theta0,theta1) {
            copClayton@paraConstr(theta0) &&
            copClayton@paraConstr(theta1) && theta1 >= theta0
        },
        ## V0 and V01
        V0 = function(n,theta) { rgamma(n, shape = 1/theta) },
        V01 = function(V0,theta0,theta1) { retstable(theta0/theta1,V0,1) },
        ## Kendall's tau
        tau = function(theta) { theta/(theta+2) },
        tauInv = function(tau) { 2*tau/(1-tau) },
        ## lower tail dependence coefficient lambda_l
        lTDC = function(theta) { 2^(-1/theta) },
        lTDCInv = function(lambda) { -1/log2(lambda) },
        ## upper tail dependence coefficient lambda_u
        uTDC = function(theta) { 0*theta },
        uTDCInv = function(lambda) {
            if(any(lambda != 0))
                stop("Any parameter for a CLayton copula gives zero upper tail dependence coefficient")
            NA * lambda
        })

stopifnot(validObject(copClayton))# ok

### ====Frank, see Nelsen (2007) p. 116, # 5====
copFrank <-
    new("ACopula", name = "Frank",
        ## generator
        psi = function(t,theta) { -log(1-(1-exp(-theta))*exp(-t))/theta },
        psiInv = function(t,theta) { -log((exp(-theta*t)-1)/(exp(-theta)-1)) },
        ## parameter interval
        paraInterval = interval("(0,Inf)"),
        ## nesting constraint
        nestConstr = function(theta0,theta1) {
            copFrank@paraConstr(theta0) &&
            copFrank@paraConstr(theta1) && theta1 >= theta0
        },
        ## V0 (algorithm of Kemp (1981)) and V01
        V0 = function(n,theta) {
            oneV0=function(theta){
              W <- runif(1)
               if(W > 1-exp(-theta)) {
                   1
               } else {
                   q <- 1-exp(-theta*runif(1))
                   if(W < q*q) { 
                     floor(1+log(W)/log(q))
                   } else if(W>q){ 
                     1 
                   } else{
                     2
                   }
               }
            }
            sapply(rep(theta,n),oneV0)
        },
        V01 = function(V0,theta0,theta1) {
            ## compute the function values of the discrete distribution up to 1-eps and maximal noValues-many values
            eps <- 1e-5
            noValues <- 500000
            Fvalues <- rep(1,noValues)
            alpha <- theta0/theta1
            c0 <- 1-exp(-theta0)
            c1 <- 1-exp(-theta1)
            ysum <- (alpha/c0)*c1       # y_1
            expsum <- log(alpha)
            logc1 <- log(c1)
            k <- 1
            while(ysum < 1-eps && k <= noValues) {
                Fvalues[k] <- ysum
                k <- k+1
                ## compute y_k
                expsum <- expsum+log(k-1-alpha)
                ysum <- ysum+exp(expsum+k*logc1-log(factorial(k)))/c0
            }
            m <- length(Fvalues)
            ## sample by finding quantiles
            n <- length(V0)
            variates <- numeric(n)
            ## warning
            numoftruncations <- 0
            for(i in 1:n) {
                uniformvariates <- runif(V0[i])
                ## sample the summands of the sums involved
                variatesfromF <- apply(outer(uniformvariates,Fvalues,">"),1,sum)+1
                numoftruncations <- numoftruncations+length(variatesfromF[variatesfromF == m+1])
                variates[i] <- sum(variatesfromF)
            }
            if(numoftruncations > 0) {
                warning("The distribution function F involved in sampling the inner distribution function for nested Frank copulas is truncated ",
                        numoftruncations," times at ",m+1,
                        " as the largest computed value of F is ",Fvalues[m])
            }
            variates
        },
        ## Kendall's tau
        ##todo: replace (1/theta)*integrate(f=function(t) t/(exp(t)-1),lower=0,upper=theta)[[1]] by debye_1(theta)
        ## tau should be: 1+(4/theta)*(D1(theta)-1), where D1 denotes the Debye function of order 1
        tau = function(theta) { sapply(theta,function(theta) { 
                                                      1+(4/theta)*((1/theta)*integrate(f=function(t) t/(exp(t)-1),
                                                      lower=0,upper=theta)[[1]]-1) })},
        tauInv = function(tau) { sapply(tau,function(tau){
                                                      r <- uniroot(function(th) copFrank@tau(th) - tau,
                                                                   interval = c(0.001,100))
                                                      ## FIXME: check for convergence
                                                      r$root
                                })          
        },
        ## lower tail dependence coefficient lambda_l
        lTDC = function(theta) { 0*theta },
        lTDCInv = function(lambda) {
            if(any(lambda != 0))
                stop("Any parameter for a Frank copula gives zero lower tail dependence coefficient")
            NA * lambda
        },
        ## upper tail dependence coefficient lambda_u
        uTDC = function(theta) { 0*theta },
        uTDCInv = function(lambda) {
            if(any(lambda != 0))
                stop("Any parameter for a Frank copula gives zero upper tail dependence coefficient")
            NA * lambda
        })

stopifnot(validObject(copFrank))# ok

### ====Gumbel, see Nelsen (2007) p. 116, # 4====
copGumbel <-
    new("ACopula", name = "Gumbel",
        ## generator
        psi = function(t,theta) { exp(-t^(1/theta)) },
        psiInv = function(t,theta) { (-log(t))^theta },
        ## parameter interval
        paraInterval = interval("[1,Inf)"),
        ## nesting constraint
        nestConstr = function(theta0,theta1) {
            copGumbel@paraConstr(theta0) &&
            copGumbel@paraConstr(theta1) && theta1 >= theta0
        },
				## V0 and V01
        V0 = function(n,theta) {
					if(theta == 1) {
				       ## Sample from S(1,1,0,1;1) with Laplace-Stieltjes transform exp(-t)
			        rep.int(1., n)
			    } else {
			        alpha <- 1/theta
			        ## Sample from S(alpha,1,(cos(alpha*pi/2))^(1/alpha),0;1)
			        ## with Laplace-Stieltjes transform exp(-t^alpha)
			        rstable1(n, alpha, beta=1,
			                 gamma = cos(alpha*pi/2)^(1/alpha))
			    }
				},
        V01 = function(V0,theta0,theta1) {
				  n <- length(V0)
			    variates <- numeric(n)
			    alpha <- theta0/theta1
			    if(alpha == 1) {
			        return(V0) # sample from S(1,1,0,V0;1) with Laplace-Stieltjes transform exp(-V0*t)
			    } else {
			        for(i in 1:n) {
			            variates[i] <- rstable1(1, alpha, beta=1,
			                                   gamma = (cos(alpha*pi/2)*V0[i])^(1/alpha))
			            ## sample from S(alpha,1,(cos(alpha*pi/2)V0)^(1/alpha),0;1) with Laplace-Stieltjes transform exp(-V0*t^alpha)
			        }
			    }
			    variates
	 			},
				## Kendall's tau
        tau = function(theta) { (theta-1)/theta },
        tauInv = function(tau) { 1/(1-tau) },
				## lower tail dependence coefficient lambda_l
        lTDC = function(theta) { 0*theta },
        lTDCInv = function(lambda) {
					if(any(lambda != 0))
              stop("Any parameter for a Gumbel copula gives zero lower tail dependence coefficient")
          NA * lambda
	 			},
        ## upper tail dependence coefficient lambda_u
        uTDC = function(theta) { 2 - 2^(1/theta) },
        uTDCInv = function(lambda) { 1/log2(2-lambda) })

stopifnot(validObject(copGumbel))# ok

## define auxiliary function for computing the mass functions for V0J and V01J
auxJoeC <- function(alpha, eps = 1e-5, noValues = 500000) {
    F <- rep(1,noValues)
    F[k <- 1] <- alpha                  # y_1
    ## MM: FIXME  vectorize -- at least in chunks;  choose() !
    while(F[k] < 1-eps && k <= noValues) {
        k <- k+1
        F[k] <- F[k-1] + choose(alpha,k)*(-1)^(k-1)
    }
    F
}

### ====Joe, see Nelsen (2007) p. 116, # 6====
copJoe <-
    new("ACopula", name = "Joe",
        ## generator
        psi = function(t,theta) { 1 - (1-exp(-t))^(1/theta) },
        psiInv = function(t,theta) { -log(1-(1-t)^theta) },
        ## parameter interval
        paraInterval = interval("[1,Inf)"),
        ## nesting constraint
        nestConstr = function(theta0,theta1) {
            copJoe@paraConstr(theta0) &&
            copJoe@paraConstr(theta1) && theta1 >= theta0
        },
        ## V0 (algorithm of Devroye (1986, p. 548)) and V01
        ## --- NB: Use auxiliary function  auxJoeC()  << FIXME give reasonable name !!
        V0 = function(n,theta) {
            ## compute the function values of the discrete distribution
            Fvalues <- auxJoeC(1/theta)
            m <- length(Fvalues)
            ## sample by finding quantiles
            uniformvariates <- runif(n)
            variates <- apply(outer(uniformvariates,Fvalues,">"),1,sum)+1
            ## warning
            numoftruncations <- length(variates[variates == m+1])
            if(numoftruncations > 0) {
                warning("The distribution function corresponding to Joe's family is truncated ",
                        numoftruncations," times at ",m+1,
                        " as the largest computed value of F is ",Fvalues[m])
            }
            variates
        },
        V01 = function(V0,theta0,theta1) {
            ## compute the function values of the discrete distribution up to 1-eps and maximal noValues-many values
            eps <- 1e-5
            noValues <- 500000
            ## compute the function values of the discrete distribution
            alpha <- theta0/theta1
            Fvalues <- auxJoeC(alpha, eps = eps, noValues = noValues)
            m <- length(Fvalues)
            ## sample by finding quantiles
            n <- length(V0)
            variates <- numeric(n)
            ## consider warnings
            n_trunc <- 0
            for(i in 1:n) {
                U.v0.i <- runif(V0[i])
                variatesfromF <- apply(outer(U.v0.i,Fvalues,">"),1,sum)+1 # sample the summands of the sums involved
                n_trunc <- n_trunc+length(variatesfromF[variatesfromF == m+1])
                variates[i] <- sum(variatesfromF)
            }
            if(n_trunc > 0) {
                warning("The distribution function F involved in sampling the inner distribution function for nested Joe copulas is truncated ",
                        n_trunc," times at ",m+1,
                        " as the largest computed value of F is", Fvalues[m])
            }
            variates
        },
        ## Kendall's tau
        tau = function(theta,noTerms=100000) { sapply(theta,
                                function(theta,noTerms=100000){
                                    k <- seq_len(noTerms)
                                    1 - 4*sum(1/(k*(theta*k+2)*(theta*(k-1)+2)))
                                })       
        },
        tauInv = function(tau) { sapply(tau,function(tau){
                        r <- uniroot(function(th) copJoe@tau(th) - tau,
                                       interval = c(1.001, 100))
                          ## FIXME: check for convergence
                          r$root
                        })
        },
        ## lower tail dependence coefficient lambda_l
        lTDC = function(theta) { 0*theta },
        lTDCInv = function(lambda) {
            if(any(lambda != 0))
                stop("Any parameter for a Joe copula gives zero lower tail dependence coefficient")
            NA * lambda
        },
        ## upper tail dependence coefficient lambda_u
        uTDC = function(theta) { 2-2^(1/theta) },
        uTDCInv = function(lambda) { log(2)/log(2-lambda) })

stopifnot(validObject(copJoe))# ok
