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
        V01 = function(V0,theta0,theta1) {
          rnbinom(length(V0),V0,(1-theta1)/(1-theta0))+V0  
        },
        ## Kendall's tau
        tau = function(theta) {
            1 - 2*((1-theta)*(1-theta)*log(1-theta)+theta)/(3*theta*theta)
        },
        tauInv = function(tau, tol = .Machine$double.eps^0.25, ...) {
            if(any(tau > 1/3))
                stop("Impossible for AMH copula to attain a Kendall's tau larger than 1/3")
            sapply(tau,function(tau) {
		r <- safeUroot(function(th) copAMH@tau(th) - tau,
			       interval = c(1e-12, 1-1e-12), Sig = +1, tol = tol, ...)
                r$root ## FIXME: check for convergence
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

## rng log(p) distribution
rlog=function(n,p){
  h=log(1-p)
  rn=numeric(n)
  for(i in 1:n){
    u2=runif(1)
    if(u2>p){
      rn[i]=1
    }else{
      u1=runif(1)
      q=1-(1-p)^u1
      if(u2<q*q){
        rn[i]=floor(1+log(u2)/log(q))
      }else if(u2>q){
        rn[i]=1
      }else{
        rn[i]=2
      }
    }
  }
  rn
}

## rng F Frank
rF=function(n,theta0,theta1){
  c=theta0/(1-exp(-theta0))
  ptheta1=1-exp(-theta1)
  alpha=theta0/theta1
  rn=numeric(n)
  for(i in 1:n){
    while(1==1){
      U=runif(1)
      X=rlog(1,ptheta1)
      if(U<=abs(choose(alpha,X))*X/alpha) break
    }
    rn[i]=X
  }
  rn
}

## Frank object
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
          p=1-exp(-theta)
          rlog(n,p)
        },
        V01 = function(V0,theta0,theta1) {
          sapply(sapply(V0,rF,theta0=theta0,theta1=theta1),sum)
        },
        ## Kendall's tau; debye_1() is from package 'gsl' :
        tau = function(theta) 1 + 4*(debye_1(theta) - 1)/theta,
	tauInv = function(tau, tol = .Machine$double.eps^0.25, ...) {
	    sapply(tau, function(tau) {
		r <- safeUroot(function(th) copFrank@tau(th) - tau,
			       interval = c(0.001,100), Sig = +1, tol=tol, ...)
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

### ====Joe, see Nelsen (2007) p. 116, # 6====

## compute Gamma constant
gammaconst=function(alpha) gamma(1-alpha)^(-1/alpha)

## define quantile function FinvJoe; gconst=gamma(1-alpha)^(-1/alpha)
FinvJoe=function(y,alpha,gconst) ifelse(y<alpha,1,floor((1-y)^(-1/alpha)*gconst))

## rng F Joe
rJoe=function(n,alpha,gconst) sapply(runif(n),FinvJoe,alpha=alpha,gconst=gconst)

## Joe object
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
        ## V0 and V01
        V0 = function(n,theta) {
          alpha=1/theta
          gconst=gammaconst(alpha)
          rJoe(n,alpha,gconst)
        },
        V01 = function(V0,theta0,theta1) {
          alpha=theta0/theta1
          gconst=gammaconst(alpha)
          sapply(sapply(V0,rJoe,alpha=alpha,gconst=gconst),sum)
        },
        ## Kendall's tau
        tau = function(theta,noTerms=100000) {
            sapply(theta,
                   function(theta,noTerms=100000) {
                       k <- seq_len(noTerms)
                       1 - 4*sum(1/(k*(theta*k+2)*(theta*(k-1)+2)))
                   })
        },
        tauInv = function(tau, tol = .Machine$double.eps^0.25, ...) {
            sapply(tau,function(tau) {
		r <- safeUroot(function(th) copJoe@tau(th) - tau,
			       interval = c(1.001, 100), Sig = +1, tol = tol, ...)
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
