##' Outer power transformation of Archimedean copulas
##'
##' @title Outer power Archimedean copulas
##' @param copbase a "base" copula, i.e. of class "acopula";
##'    must be one of the 5 five predefined families
##' @param thetabase the (univariate) parameter 'theta' for the base copula
##' @return a new "acopula" object; the outer power copula
##' @author Marius Hofert
opower <- function(copbase, thetabase) {
    ## create object with name in here so it's in the environment and we can access it
    cOP <- new("acopula", name = paste("opower", copbase@name, sep=":"),
               ## generator
               psi = function(t,theta) { copbase@psi(t^(1/theta), thetabase) },
               psiInv = function(t,theta) { copbase@psiInv(t, thetabase)^theta },
               ## parameter interval
               paraInterval = interval("[1,Inf)"),
               ## nesting constraint
               nestConstr = function(theta0,theta1) {
                   copbase@paraConstr(theta0) &&
                   copbase@paraConstr(theta1) && theta1 >= theta0
               },
               ## absolute value of generator derivatives
               psiDabs = function(t, theta, degree=1, n.MC=0, log=FALSE){
                   if(theta == 1) return(copbase@psiDabs(t, theta, degree=degree, n.MC=n.MC, log=log)) # copbase case
                   if(n.MC > 0){
                       psiDabsMC(t, family=cOP, theta=theta, degree=degree, n.MC=n.MC, log=log)
                   }else{
                       ## FIXME: not optimal yet, inner sum could get a *real* log
                       j <- 1:degree
                       facs <- (-1)^j*unlist(lapply(j/theta, function(z) prod(z-(0:(degree-1)))))
                       t.beta <- t^(1/theta)
                       l.t.beta <- (1/theta)*log(t)
                       fun <- function(x,y) copbase@psiDabs(y, thetabase, degree=x, n.MC=n.MC, log=TRUE)
                       l.psiDabs.t.beta <- outer(j, t.beta, FUN=fun)
                       one.j <- function(j.){
                           k <- j.:degree
                           colSums(exp(outer(k-1, l.t.beta)-lfactorial(j.)-lfactorial(k-j.)+l.psiDabs.t.beta[k,]))
                       }
                       mat <- do.call(rbind, lapply(j, one.j)) # degree x length(t) matrix
                       facs %*% mat/t^(degree-1/theta)
                   }
               },
               ## derivatives of the generator inverse
               psiInvD1abs = function(t, theta, log=FALSE){
                   if(theta == 1) return(copbase@psiInvD1abs(t, theta, log=log)) # copbase case
                   if(log){
                       log(theta)+(theta-1)*log(copbase@psiInv(t,thetabase))+
                           copbase@psiInvD1abs(t, thetabase,log=TRUE)
                   }else{
                       theta*copbase@psiInv(t,thetabase)^(theta-1)*
                           copbase@psiInvD1abs(t, thetabase,log=FALSE)
                   }
               },
               ## density
               dacopula = function(u, theta, n.MC=0, log=FALSE){
                   if(theta == 1) return(copbase@dacopula(u, theta, n.MC=n.MC, log=log)) # copbase case
                   if(!is.matrix(u)) u <- rbind(u)
                   if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
                   ## f() := NaN outside and on the boundary of the unit hypercube
                   res <- rep.int(NaN, n <- nrow(u))
                   n01 <- apply(u,1,function(x) all(0 < x, x < 1)) # indices for which density has to be evaluated
                   if(!any(n01)) return(res)
                   ## auxiliary results
                   u. <- u[n01,, drop=FALSE]
                   psiI <- rowSums(cOP@psiInv(u.,theta))
                   res[n01] <- cOP@psiDabs(psiI, theta, degree=d, n.MC=n.MC, log=TRUE)+
                       rowSums(cOP@psiInvD1abs(u., theta, log=TRUE))
                   if(log) res else exp(res)
               },
               ## V0 and V01
               V0 = function(n,theta) {
                   if(theta == 1) {
                       ## Sample from S(1,1,0,1;1)
                       ## with Laplace-Stieltjes transform exp(-t)
                       rep.int(1., n)
                   } else {
                       V0base <- copbase@V0(n,thetabase) # draw from the base generator
                       alpha <- 1/theta
                       ## Sample from S(alpha,1,(cos(alpha*pi/2))^(1/alpha),0;1)
                       ## with Laplace-Stieltjes transform exp(-t^alpha)
                       S <- rstable1(n, alpha, beta=1,
                                     gamma = (cos(alpha*pi/2))^(1/alpha))
                       S*V0base^theta
                   }
               },
               dV0 = function(x, theta, log=FALSE){
                   stop("not implemented; it's the density of SV^theta, where V ~ F with LS[F] = copbase@psi and S ~ S(1/theta, 1, cos^theta(pi/(2*theta)), I_{theta==1}; 1)")
               },
               V01 = function(V0,theta0,theta1) {
                   alpha <- theta0/theta1
                   if(alpha == 1) {
                       ## Sample from S(1,1,0,V0;1)
                       ## with Laplace-Stieltjes transform exp(-V0*t)
                       V0
                   } else {
                       rstable1(length(V0), alpha, beta=1,
                                gamma = (cos(alpha*pi/2)*V0)^(1/alpha))
                       ## Sample from S(alpha,1,(cos(alpha*pi/2)V0)^(1/alpha),0;1)
                       ## with Laplace-Stieltjes transform exp(-V0*t^alpha)
                   }
               },
               dV01 = function(x, V0, theta0, theta1, log=FALSE){
                   copGumbel@dV01(x, V0, theta0, theta1, log=log)
               },
               ## Kendall's tau
               tau = function(theta) {
                   1-(1-copbase@tau(thetabase))/theta
               },
               tauInv = function(tau) {
                   taubase <- copbase@tau(thetabase)
                   if(tau >= taubase) (1-taubase)/(1-tau)
                   else {
                       stop("The provided tau has to be >= taubase")
                       NA * tau
                   }
               },
               ## lower tail dependence coefficient lambda_l
               lambdaL = function(theta) {
                   if(copbase@name=="Clayton") 2^(-1/(thetabase*theta)) else 0*theta
               },
               lambdaLInv = function(lambda) {
                   if(copbase@name=="Clayton") {
                       if(lambda >= 2^(-1/thetabase)) -1/(thetabase*log2(lambda))
                       else {
                           stop("The provided lambda has to be >= 2^(-1/thetabase)")
                           NA * lambda
                       }
                   } else {
                       if(any(lambda != 0))
                           stop("Any parameter for this outer power Archimedean copula gives lambdaL = 0")
                       NA * lambda
                   }
               },
               ## upper tail dependence coefficient lambda_u
               lambdaU = function(theta) {
                   2 - 2^(1/if(copbase@name %in% c("Gumbel", "Joe")) thetabase*theta else theta)
               },
               lambdaUInv = function(lambda) {
                   if(copbase@name %in% c("Gumbel", "Joe")) {
                       if(lambda >= 2-2^(1/thetabase)) 1/(thetabase*log2(2-lambda))
                       else {
                           stop("The provided lambda has to be >= 2-2^(1/thetabase)")
                           NA * lambda
                       }
                   } else 1/log2(2-lambda)
               }
               )
    cOP
}
