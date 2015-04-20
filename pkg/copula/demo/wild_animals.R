require(lattice)

### 1) Example from Marius Hofert and Frederik Vrins ###########################

## Lambda and its inverse
Lambda <- function(t) pmax(log(t), 0)
LambdaInv <- function(t) ifelse(t == 0, 0, exp(t))

## M and its inverse (for M_i)
M1 <- function(t) 0.01*t
M1Inv <- function(t) t/0.01
M2 <- function(t){
    f <- function(t.){
	ct <- ceiling(t.)
	if(ct %% 2) t.-(ct-1)/2 else ct/2
    }
    unlist(lapply(t, f))
}
M2Inv <- function(t){
    f <- function(t.){
        if(t. == 0) return(0)
        t. + ceiling(t) - 1
    }
    unlist(lapply(t, f))
}

## S and its inverse (for S_i)
S1 <- function(t, H) exp(-(M1(t) + Lambda(t)*(1-exp(-H))))
S1Inv <- function(u, H, upper=1e6) unlist(lapply(u, function(u.)
    uniroot(function(x) S1(x,H=H)-u., interval=c(0, upper))$root))
S2 <- function(t, H) exp(-(M2(t) + Lambda(t)*(1-exp(-H))))
S2Inv <- function(u, H, upper=1e6) unlist(lapply(u, function(u.)
    uniroot(function(x) S2(x,H=H)-u., interval=c(0, upper))$root))

## wrappers for p_1 and p_2 and their inverses
p <- function(t, k, H, I, p1, p2){
    if((lI <- length(I)) == 0){
        stop("error in p")
    }else if(lI==1){
        if(I==1) p1(t, k=k, H=H) else p2(t, k=k, H=H)
    }else{ # lI == 2
        c(p1(t, k=k, H=H), p2(t, k=k, H=H))
    }
}
pInv <- function(u, k, H, I, p1Inv, p2Inv){
    if((lI <- length(I)) == 0){
        stop("error in pInv")
    }else if(lI==1){
        if(I==1) p1Inv(u, k=k, H=H) else p2Inv(u, k=k, H=H)
    }else{ # lI == 2
        c(p1Inv(u[1], k=k, H=H), p2Inv(u[2], k=k, H=H))
    }
}

## Define C
C <- function(u, H, Lambda, S1Inv, S2Inv){
    if(u[1]==0 & u[2]==0) 0 else
    u[1]*u[2]*exp(expm1(-H)^2*Lambda(min(S1Inv(u[1],H), S2Inv(u[2],H))))
}

## Compute the singular component (given u1, find u2 on the singular component)
## u2=S2(S1^-(u1))
s.comp <- function(u1, H, S1Inv, S2){
    if(u1 == 0) 0
    else S2(S1Inv(u1,H),H)
}

## p and its inverse (for p_i(t_k-) and p_i(t_k), i = 1,2)
p1 <- function(t, k, H) exp(-M1(t)-H*k)
p1Inv <- function(u, k, H, upper=1e6) unlist(lapply(u, function(u.)
    uniroot(function(x) p1(x,k=k,H=H)-u., interval=c(0, upper))$root))
p2 <- function(t, k, H) exp(-M2(t)-H*k)
p2Inv <- function(u, k, H, upper=1e6) unlist(lapply(u, function(u.)
    uniroot(function(x) p2(x,k=k,H=H)-u., interval=c(0, upper))$root))

## Generate one bivariate random vector from C
rC1 <- function(H, LambdaInv, S1, S2, p1, p2, p1Inv, p2Inv){
    d <- 2 # dim = 2
    ## (1)
    U <- runif(d) # for determining the default times of all components
    ## (2)
    t_h <- 0 # t_{h,0}; initial value for the occurrence of the homogeneous Poisson process with unit intensity
    t <- 0 # t_0; initial value for the occurrence of the jump process
    k <- 1 # indices for the sets I
    I_ <- list(1:d) # I_k; indices for which tau has to be determined
    tau <- rep(Inf,d) # in the beginning, set all default times to Inf (= "no default")
    ## (3)
    repeat{
        ## (4)
        t_h[k] <- rexp(1) + if(k == 1) 0 else t_h[k-1] # k-th occurrence of a homogeneous Poisson process with unit intensity
        t[k] <- LambdaInv(t_h[k]) # k-th occurrence of a non-homogeneous Poisson process with integrated rate function Lambda
        ## (5)
        pvec <- p(t[k], k=k, H=H, I=c(1,2), p1=p1, p2=p2) # c(p_1(t_k), p_2(t_k))
        I. <- (1:d)[U >= pvec] # determine all i in I
        I <- intersect(I., I_[[k]]) # determine I
        ## (6)--(10)
        if(length(I) > 0){
            default.at.t <- U[I] <= p(t[k], k=k-1, H=H, I=I, p1=p1, p2=p2)
            tau[I[default.at.t]] <- t[k]
            Ic <- I[!default.at.t] # I complement
            if(length(Ic) > 0) tau[Ic] <- pInv(U[Ic], k=k-1, H=H, I=Ic,
                                               p1Inv=p1Inv, p2Inv=p2Inv)
        }
        ## (11)
        I_[[k+1]] <- setdiff(I_[[k]],I) # define I_{k=1} = I_k \ I
        ## (12)
        if(length(I_[[k+1]]) == 0) break else k <- k+1
    }
    ## (14)
    c(S1(tau[1], H=H), S2(tau[2], H=H))
}

## Draw n vectors of random variates from C
rC <- function(n, H, LambdaInv, S1, S2, p1, p2, p1Inv, p2Inv){
    mat <- t(sapply(rep(H, n), rC1, LambdaInv=LambdaInv, S1=S1, S2=S2,
                    p1=p1, p2=p2, p1Inv=p1Inv, p2Inv=p2Inv))
    row.names(mat) <- NULL
    mat
}

## Generate copula data
n <- 2e4
H <- 10
set.seed(271)
U <- rC(n, H=H, LambdaInv=LambdaInv, S1=S1, S2=S2, p1=p1, p2=p2, p1Inv=p1Inv, p2Inv=p2Inv)

## Check margins of U
hist(U[,1], probability=TRUE, main="Histogram of the first component",
     xlab=expression(italic(U[1])))
hist(U[,2], probability=TRUE, main="Histogram of the second component",
     xlab=expression(italic(U[2])))

## Plot U (copula sample)
plot(U, cex=0.2, xlab=expression(italic(U[1])%~%~"U[0,1]"),
                 ylab=expression(italic(U[2])%~%~"U[0,1]"))

## Wireframe plot to incorporate singular component
wf.plot <- function(grid, val.grid, s.comp, val.s.comp, Lambda, S1Inv, S2Inv){
    wireframe(val.grid ~ grid[,1]*grid[,2], xlim=c(0,1), ylim=c(0,1), zlim=c(0,1),
              aspect=1, scales=list(arrows=FALSE, col=1), # remove arrows
              par.settings=list(axis.line=list(col="transparent"), # remove global box
              clip=list(panel="off")), pts=cbind(s.comp, val.s.comp), # add singular component
              panel.3d.wireframe = function(x,y,z,xlim,ylim,zlim,xlim.scaled,
              ylim.scaled,zlim.scaled,pts,...){
                  panel.3dwire(x=x,y=y,z=z,xlim=xlim,ylim=ylim,zlim=zlim,
                               xlim.scaled=xlim.scaled,ylim.scaled=ylim.scaled,
                               zlim.scaled=zlim.scaled,...)
                  xx <- xlim.scaled[1]+diff(xlim.scaled)*(pts[,1]-xlim[1])/diff(xlim)
                  yy <- ylim.scaled[1]+diff(ylim.scaled)*(pts[,2]-ylim[1])/diff(ylim)
                  zz <- zlim.scaled[1]+diff(zlim.scaled)*(pts[,3]-zlim[1])/diff(zlim)
                  panel.3dscatter(x=xx,y=yy,z=zz,xlim=xlim,ylim=ylim,zlim=zlim,
                                  xlim.scaled=xlim.scaled,ylim.scaled=ylim.scaled,
                                  zlim.scaled=zlim.scaled,type="l",col=1,...)
              },
              xlab=expression(italic(u[1])),
              ylab=expression(italic(u[2])),
              zlab=list(expression(italic(C(u[1],u[2]))), rot=90))
}

## Copula plot with singular component
u <- seq(0, 1, length.out=20) # grid points per dimension
grid <- expand.grid(u1=u, u2=u) # grid
val.grid <- apply(grid, 1, C, H=H, Lambda=Lambda, S1Inv=S1Inv, S2Inv=S2Inv) # copula values on grid
s.comp <- cbind(u, sapply(u, s.comp, H=H, S1Inv=S1Inv, S2=S2)) # pairs (u1, u2) on singular component
val.s.comp <- apply(s.comp, 1, C, H=H, Lambda=Lambda, S1Inv=S1Inv, S2Inv=S2Inv) # corresponding z-values
wf.plot(grid=grid, val.grid=val.grid, s.comp=s.comp, val.s.comp=val.s.comp,
        Lambda=Lambda, S1Inv=S1Inv, S2Inv=S2Inv)


### 2) A nice example from Wolfgang Trutschnig and Manuela Schreyer ############

## See Trutschnig, Fernandez Sanchez (2014, "Copulas with continuous, strictly
## increasing singular conditional distribution functions") for more details

## Roughly, one defines an Iterated Function System whose attractor is the word
## 'Copula' and starts the chaos game

## Define the Iterated Function System
n <- 23
IFS <- list(function(x) c(3*x[1]/n,      x[2]/4),
            function(x) c(-(x[2]-1)/n,   x[1]/2+1/4),
            function(x) c(3*x[1]/n,      x[2]/4+3/4),
            function(x) c((3*x[1]+4)/n,  x[2]/4),
            function(x) c(-(x[2]-5)/n,   x[1]/2+1/4),
            function(x) c((3*x[1]+4)/n,  x[2]/4+3/4),
            function(x) c(-(x[2]-7)/n,   x[1]/2+1/4),
            function(x) c(-x[2]/n+9/n,   3*x[1]/4),
            function(x) c((3*x[1]+8)/n,  x[2]/4+3/4),
            function(x) c(x[1]/n+10/n,   x[2]/8+1/2+1/8),
            function(x) c(2*x[1]/n+9/n,  x[2]/4+1/4+1/8),
            function(x) c(-x[2]/n+13/n,  (3*x[1]+1)/4),
            function(x) c((3*x[1]+12)/n, x[2]/4),
            function(x) c((3*x[1]+12)/n, x[2]/4),
            function(x) c(-x[2]/n+15/n,  (3*x[1]+1)/4),
            function(x) c((3*x[1]+16)/n, x[2]/4),
            function(x) c(-(x[2]-21)/n,  3*x[1]/4),
            function(x) c((3*x[1]+20)/n, x[2]/4+3/4),
            function(x) c((x[1]+21)/n,   x[2]/4+1/4+1/8),
            function(x) c(-(x[2]-23)/n,  3*x[1]/4))

## Run chaos game B times
B <- 20 # replications
n.steps <- 20000 # number of steps
AA <- vector("list", length=B)
set.seed(271)
for(i in 1:B) {
    ind <- sample(length(IFS), size=n.steps, replace=TRUE) # (randomly) 'bootstrap' functions
    res <- matrix(0, nrow=n.steps+1, ncol=2) # result matrix (for each i)
    pt <- c(0, 0) # initial point
    for(r in seq_len(n.steps)) {
        res[r+1,] <- IFS[[ind[r]]](pt) # evaluate randomly chosen functions at pt
        pt <- res[r+1,] # redefine point
    }
    AA[[i]] <- res # keep this matrix
}
A <- do.call("rbind", AA) # rbind (n.steps+1, 2)-matrices
stopifnot(nrow(A) == B*(n.steps+1)) # sanity check

## Rotate (=> X)
phi <- -pi/4
X <- cbind(cos(phi)*A[,1]-sin(phi)*A[,2]/3,
           sin(phi)*A[,1]+cos(phi)*A[,2]/3)

## Transform with marginal edfs (=> U)
U <- apply(X, 2, function(x) ecdf(x)(x))

## Check margins of U (here 'perfect')
hist(U[,1], probability=TRUE, main="Histogram of the first component",
     xlab=expression(italic(U[1])))
hist(U[,2], probability=TRUE, main="Histogram of the second component",
     xlab=expression(italic(U[2])))

## Plot U (copula sample)
plot(U, pch=".", xlab=expression(italic(U[1])%~%~"U[0,1]"),
                 ylab=expression(italic(U[2])%~%~"U[0,1]"))
