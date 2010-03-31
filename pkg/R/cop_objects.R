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

## rng Log(p) distribution
rlog=function(n,p){
  vec=numeric(n)
  u=runif(n)
  l1=(u>p)
  i1=(1:n)[l1]
  vec[i1]=1
  l2=!l1
  i2=(1:n)[l2]
  q=1-(1-p)^runif(n)
  l3=(u[i2]<q[i2]*q[i2])
  i3=i2[l3]
  vec[i3]=floor(1+log(u[i3])/log(q[i3]))
  l4=(u[i2]>q[i2])
  i4=i2[l4]
  vec[i4]=rep(1,length(i4))
  l5=!(l3|l4)
  i5=i2[l5]
  vec[i5]=rep(2,length(i5))
  vec 
}

## rejection for F for Frank's family
rejFFrank=function(p,alpha,theta0le1){
  if(theta0le1){
    while(1==1){
      u=runif(1)
      x=rlog(p)
      if(u<=1/((x-alpha)*beta(x,1-alpha))) break
    }
  }else{
    while(1==1){
      u=runif(1)
      x=rFJoe(1,alpha)
      if(u<=p^(x-1)) break
    }
  }
  x
}

## sample F for Frank 
rFFrank=function(n,theta0,theta1){ sapply(rep(1-exp(-theta1),n),rejFFrank,alpha=theta0/theta1,theta0le1=ifelse(theta0<=1,TRUE,FALSE)) }

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
        V0 = function(n,theta) { rlog(n,1-exp(-theta)) },
        V01 = function(V0,theta0,theta1) {
          sapply(lapply(V0,rFFrank,theta0=theta0,theta1=theta1),sum) ##FIXME: hoe to approximate when V0 large? not solved yet
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
            alpha <- theta0/theta1
            if(alpha == 1) {
                V0 # sample from S(1,1,0,V0;1) with Laplace-Stieltjes transform exp(-V0*t)
            } else {
                rstable1(length(V0), alpha, beta=1,
                         gamma = (cos(alpha*pi/2)*V0)^(1/alpha))
                ## sample from S(alpha,1,(cos(alpha*pi/2)V0)^(1/alpha),0;1) with Laplace-Stieltjes transform exp(-V0*t^alpha)
            }
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

rFJoe=function(n,alpha){
  vec=numeric(n)
  u=runif(n)
  l1=(u<=alpha)
  i1=(1:n)[l1]
  vec[i1]=1
  l2=!l1
  i2=(1:n)[l2]
  Ginv=((1-u[i2])*gamma(1-alpha))^(-1/alpha)
  floorGinv=floor(Ginv)
  l3=(1-1/(floorGinv*beta(floorGinv,1-alpha))<u[i2])
  i3=i2[l3]
  vec[i3]=ceiling(Ginv[l3])
  l4=!l3
  i4=i2[l4]
  vec[i4]=floorGinv[l4]
  vec
}

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
        V0 = function(n,theta) { rFJoe(n,1/theta) },
        V01 = function(V0,theta0,theta1,approx=100000) {
          alpha=theta0/theta1
          ia=(V0>approx)
          ie=!ia
          V01=numeric(length(V0))
          V01[ia]=V0[ia]^(1/alpha)*rstable1(sum(ia),alpha,1,(cos(alpha*pi/2))^(1/alpha),ifelse(alpha==1,1,0))
          V01[ie]=sapply(lapply(V0[ie],rFJoe,alpha=alpha),sum) 
          V01
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
