##' Transforms vectors of random variates following the given (nested) Archimedean
##' copula (with specified parameters) to U[0,1]^d vectors of random variates
##' via Rosenblatt's transformation
##'
##' @title Rosenblatt transformation
##' @param u data matrix in [0,1]^(n, d) ((pseudo-/copula-)observations if
##'        inverse==TRUE and U[0,1] observations if inverse==FALSE)
##' @param cop object of class Copula
##' @param j.ind single index in {2,..,d} for which C(u_j | u_1,..,u_{j-1}) is
##'        computed and returned or NULL for which C(u_j | u_1,..,u_{j-1}) is
##'        computed for all j in {2,...,d} (and in which case the non-transformed
##'        first column is also returned)
##' @param n.MC parameter n.MC for evaluating the derivatives via Monte Carlo;
##'        if 0, the available (theoretical) formula is used.
##' @param inverse logical indicating whether the inverse of rtrafo is computed
##'        (this is known as 'conditional distribution method' for sampling)
##' @param log logical indicating whether the log-transform is computed
##' @return An (n,d) matrix U of supposedly U[0,1]^d realizations (if inverse==FALSE)
##'         or copula samples (if inverse=TRUE) (if j.ind==NULL) or
##'         C(u_j | u_1,..,u_{j-1}) (if j.ind in {2,..,d}) [or the log of
##'         the result if log=TRUE]
##' @author Marius Hofert and Martin Maechler
rtrafo <- function(u, cop, j.ind=NULL, n.MC=0, inverse=FALSE, log=FALSE)
{
    ## checks
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    stopifnot(0 <= u, u <= 1, is(cop, "Copula"), n.MC >= 0)
    n <- nrow(u)
    d <- ncol(u)
    stopifnot(d >= 2)
    is.null.j.ind <- is.null(j.ind)
    stopifnot(is.null.j.ind || (length(j.ind)==1 && 2 <= j.ind && j.ind <= d)) # either NULL or of length 1
    jj     <- if(is.null.j.ind) 2:d else j.ind
    jj.res <- if(is.null.j.ind) 1:d else j.ind

    ## Distinguish the families (as they are quite different)

    ## (nested) Archimedean case ###############################################

    if((NAC <- is(cop, "outer_nacopula")) || is(cop, "archmCopula")) {
	if(NAC) {
	    if(length(cop@childCops))
		stop("Currently, only Archimedean copulas are supported")
            ## outer_nacopula but with no children => an AC => continue
	    cop <- cop@copula # class(cop) = "acopula"
	    th <- cop@theta
	} else { # "archmCopula"
	    th <- cop@parameters
	    cop <- getAcop(cop) # class(cop) = "acopula" but without parameter or dim
	}
	stopifnot(cop@paraConstr(th))

        ## Compute inverse
        if(inverse) {
            U <- u # u's are supposedly U[0,1]^d
            max.col <- if(is.null(j.ind)) d else j.ind # maximal column index up to which the conditional copula inverses need to be computed
            if(cop@name=="Clayton") { # Clayton case (explicit)
                sum. <- U[,1]^(-th)
                for(j in 2:max.col) {
                    U[,j] <- log1p((1-j+1+sum.)*(u[,j]^(-1/(j-1+1/th)) - 1))/(-th)
                    eUj <- exp(U[,j])
                    sum. <- sum. + eUj^(-th)
                    if(!log) U[,j] <- eUj
                }
            } else { # non-Clayton / general case
                if(n.MC) stop("The inverse of Rosenblatt's transformation with uniroot() requires n.MC=0")
                ## TODO: After the acopula and archmCopula classes are better merged, the
                ##       tedious conversion to acopula (above) and then again to archmCopula (below)
                ##       is not required anymore.
                arCop <- archmCopula(cop@name, param=th)
                ## f(x) := cond_cop_given_j_minus_1(x) -  u[i,j]  =  rtrafo(...) - u[i,j]
                fx.ij <- function(x) rtrafo(c(U.j1[i, , drop=FALSE], x), cop=arCop, j.ind=j) - uj[i]
                for(j in 2:max.col) {
                    ## precompute j-stuff so the uniroot() in for(i ..) is as fast as possible:
                    U.j1 <- U[, seq_len(j-1L), drop=FALSE]
                    uj <- u[,j]
                    for(i in 1:n)
                        U[i,j] <- uniroot(fx.ij, interval=c(0,1))$root
                }
                if(log) U <- log(U)
            }
            return(U[,jj.res])
        }

        ## Compute conditional probabilities
	psiI <- cop@iPsi(u, theta=th)	   # (n,d) matrix of psi^{-1}(u)
	psiI. <- t(apply(psiI, 1, cumsum)) # corresponding (n,d) matrix of row sums
	if(n.MC==0){
	    ## Note: C(u_j | u_1,...,u_{j-1}) = \psi^{(j-1)}(\sum_{k=1}^j \psi^{-1}(u_k)) /
            ##                                  \psi^{(j-1)}(\sum_{k=1}^{j-1} \psi^{-1}(u_k))
	    C.j <- function(j) {
                ## computes C(u_j | u_1,...,u_{j-1}) with the same idea as cacopula()
                ## but faster
		logD <- cop@absdPsi(as.vector(psiI.[,c(j,j-1)]), theta=th,
				    degree=j-1, n.MC=0, log=TRUE)
		res <- logD[1:n]-logD[(n+1):(2*n)]
		if(log) res else exp(res)
	    }
	} else { ## n.MC > 0
	    ## draw random variates
	    V <- cop@V0(n.MC, th)
	    C.j <- function(j) {
                ## computes C(u_j | u_1,...,u_{j-1}) with the same idea as
                ## default method of absdPsiMC (only difference: draw V's only once)
		arg <- c(psiI.[,j], psiI.[,j-1])
		iInf <- is.infinite(arg)
		logD <- numeric(2*n)
		logD[iInf] <- -Inf
		if(any(!iInf)) logD[!iInf] <- lsum(-V %*% t(arg[!iInf]) +
						   (j-1) * log(V) - log(n.MC))
		res <- logD[1:n]-logD[(n+1):(2*n)]
		if(log) res else exp(res)
	    }
	}
        trafo <- vapply(jj, C.j, numeric(n))
        if(is.null.j.ind) cbind(trafo, if(log) log(u[,1]) else u[,1]) else trafo

    } else if(is(cop, "normalCopula")) {

    ## Gauss copula (see, e.g., Cambou, Hofert, Lemieux) #######################

        P <- getSigma(cop)
        stopifnot(dim(P) == c(d,d)) # defensive programming

        ## compute inverse
        if(inverse) {
            U <- u # consider u as U[0,1]^d
            max.col <- if(is.null.j.ind) d else j.ind
            x <- qnorm(u[,1:max.col, drop=FALSE]) # will be updated with previously transformed U's
            for(j in 2:max.col) {
                P. <- P[j,1:(j-1), drop=FALSE] %*% solve(P[1:(j-1),1:(j-1), drop=FALSE]) # (1,j-1) %*% (j-1,j-1) = (1,j-1)
                mu.cond <- as.numeric(P. %*% t(x[,1:(j-1), drop=FALSE])) # (1,j-1) %*% (j-1,n) = (1,n) = n
                P.cond <- P[j,j] - P. %*% P[1:(j-1),j, drop=FALSE] # (1,1) - (1,j-1) %*% (j-1,1) = (1,1)
                U[,j] <- pnorm(qnorm(u[,j], mean=mu.cond, sd=sqrt(P.cond)), log.p=log)
                x[,j] <- qnorm(if(log) exp(U[,j]) else U[,j]) # based on previously transformed U's
            }
            if(log) U[,1] <- log(U[,1]) # adjust the first column for 'log'
            return(U[,jj.res])
        }

        ## compute conditional probabilities (more efficient due to vapply())
        max.col.ran <- if(is.null.j.ind) jj.res else 1:j.ind
        x <- qnorm(u[, max.col.ran, drop=FALSE])
        C.j <- function(j) {
            P. <- P[j,1:(j-1), drop=FALSE] %*% solve(P[1:(j-1),1:(j-1), drop=FALSE]) # (1,j-1) %*% (j-1,j-1) = (1,j-1)
            mu.cond <- as.numeric(P. %*% t(x[,1:(j-1), drop=FALSE])) # (1,j-1) %*% (j-1,n) = (1,n) = n
            P.cond <- P[j,j] - P. %*% P[1:(j-1),j, drop=FALSE] # (1,1) - (1,j-1) %*% (j-1,1) = (1,1)
            pnorm(x[,j], mean=mu.cond, sd=sqrt(P.cond), log.p=log)
        }
        trafo <- vapply(jj, C.j, numeric(n))
        if(is.null.j.ind) cbind(trafo, if(log) log(u[,1]) else u[,1]) else trafo

    } else if(is(cop, "tCopula")) {

    ## t Copula (see, e.g., Cambou, Hofert, Lemieux) ###########################

	P <- getSigma(cop)
        stopifnot(dim(P) == c(d,d)) # defensive programming
	nu <- getdf(cop)
        n <- nrow(u)

        ## compute inverse
        if(inverse) {
            U <- u # consider u as U[0,1]^d
            max.col <- if(is.null(j.ind)) d else j.ind
            x <- qt(u[,1:max.col, drop=FALSE], df=nu) # will be updated with previously transformed U's
            for(j in 2:max.col) {
                P1.inv <- solve(P[1:(j-1),1:(j-1), drop=FALSE])
                x1 <- x[,1:(j-1), drop=FALSE]
                g  <- vapply(1:n, function(i) x1[i, ,drop=FALSE] %*% P1.inv %*%
                             t(x1[i, ,drop=FALSE]), numeric(1))
                P.inv <- solve(P[1:j, 1:j, drop=FALSE])
                s1 <- sqrt((nu+j-1)/(nu+g))
                s2 <- (x1 %*% P.inv[1:(j-1),j, drop=FALSE]) / sqrt(P.inv[j,j])
                U[,j] <- pt((qt(u[,j], df=nu+j-1)/s1-s2) / sqrt(P.inv[j,j]),
                            df=nu, log.p=log)
                x[,j] <- qt(if(log) exp(U[,j]) else U[,j], df=nu) # based on previously transformed U's
            }
            if(log) U[,1] <- log(U[,1]) # adjust the first column for 'log'
            return(U[,jj.res])
        }

        ## compute conditional probabilities (more efficient due to vapply())
        max.col.ran <- if(is.null.j.ind) jj.res else 1:j.ind
        x <- qt(u[, max.col.ran, drop=FALSE], df=nu)
        C.j <- function(j) {
            P1.inv <- solve(P[1:(j-1),1:(j-1), drop=FALSE])
            x1 <- x[,1:(j-1), drop=FALSE]
            g  <- vapply(1:n, function(i) x1[i, ,drop=FALSE] %*% P1.inv %*%
                         t(x1[i, ,drop=FALSE]), numeric(1))
            P.inv <- solve(P[1:j, 1:j, drop=FALSE])
            s1 <- sqrt((nu+j-1)/(nu+g))
            s2 <- (x1 %*% P.inv[1:(j-1),j, drop=FALSE]) / sqrt(P.inv[j,j])
            lres <- pt(s1 * ( sqrt(P.inv[j,j]) * x[,j, drop=FALSE] + s2),
                       df=nu+j-1, log.p=TRUE)
            if(log) lres else exp(lres)
        }
        trafo <- vapply(jj, C.j, numeric(n))
        if(is.null.j.ind) cbind(trafo, if(log) log(u[,1]) else u[,1]) else trafo

    } else {
	stop("Not yet implemented for copula class ", class(cop))
    }
}
