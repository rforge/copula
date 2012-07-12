## test of fitCopula
######################################

## one replicate
do1 <- function(cop,n) {
  x <- rcopula(cop,n)
  u <- pobs(x)
  f.itau <- fitCopula(cop, u, method="itau")
  f.irho <- fitCopula(cop, u, method="irho")
  f.mpl <- fitCopula(cop, u, method="mpl")
  f.ml <- fitCopula(cop, x, method="ml")
  c(itau.est = f.itau@estimate, irho.est = f.irho@estimate,
    mpl.est  = f.mpl @estimate, ml.est   = f.ml  @estimate,
    itau.se = sqrt(f.itau@var.est), irho.se = sqrt(f.irho@var.est),
    mpl.se  = sqrt(f.mpl @var.est), ml.se   = sqrt(f.ml  @var.est))
}

## N replicates with mean of estimates and se, and sd of estimates
testCop <- function(cop, n, N) {
  sim <- replicate(N, do1 (cop,n))
  m <- rowMeans(sim)
  cbind(bias = m[1:4] - cop@parameters, ## bias
        dSE  = m[5:8] - apply(sim[1:4,], 1, sd)) ## mean of se minus sd of estimates
}


## test for one bivariate one-parameter copula family
run1test <- function(cop, taus, n, N) {
  theta <- structure(calibKendallsTau(cop,taus),
		     names = paste("tau", format(taus), sep="="))
  setPar <- function(cop, value) { cop@parameters <- value ; cop }
  sapply(lapply(theta, setPar, cop=cop),
	 testCop, n=n, N=N,
	 simplify = "array")
}

rr <- run1test(normalCopula(dim=2), taus=seq(0.2,0.8,by=0.2), n=16, N= 8)
