#### Experiments with retstableR()

library(nacopula)

### investigate the  m(V0) = "optimal m" : ----------------------

V0 <- seq(0, 40, by = 1/256)
nm <- length(m <- sapply(V0, nacopula:::m.opt.retst))
iJmp <- which(diff(m) != 0)
stopifnot( all.equal(m[c(iJmp,nm)], 1:40) )
iKept <- c(1, rbind(iJmp, iJmp+1), nm)

V0. <- V0[iKept]
 m. <-  m[iKept]

V0[iJmp] - m[iJmp]     ## 0.3828125 0.4296875 0.4492188 0.4609375 .... converging to 0.5
V0[iJmp+1] - m[iJmp+1] ## .... converging to -0.5


### ------------ How many times is rstable1(1, ...) called ? ---------------------
### --> we *really* should call  rstable1(N, ....)
### and then use these values  where  N = N(V0):

saveFile <- "retstable_Nstat.rda"

## This needs two and a half hour -- elapsed time !!
if(file.exists(saveFile)) {
    load      (saveFile)
} else {
  .e <- new.env()
  ## our tracer function is just a "counter":
  rstTR <- function() .e$N <- if(is.null(.e$N)) 0L else .e$N + 1L
  trace(rstable1, rstTR, print=FALSE)

  ## For now, we use the R (reliable / correct) version:

  nV <- length(V0s <- c(1:4, 5*(1:10), 10*(6:20)))
  nAlph <- length(alphas <- (1:9)/10)
  nSim <- 1024

  Nstat <- array(dim = c(nSim, nAlph, nV))
  for(iV0 in 1:nV) {
      V0 <- V0s[iV0]
      cat(sprintf("V0 = %5d :", V0))
      for(ia in 1:nAlph) {
          alpha <- alphas[ia]
          Nstat[, ia, iV0] <- replicate(nSim,
                                    { .e$N <- 0L # start counting ...
                                      RR <- retstableR(alphas[ia], V0 = V0)
                                      .e$N }
                                        )
          cat(".")
      }
      cat("\n")
  }

  save(V0s,alphas, nSim, Nstat, file = saveFile)
}

## Visualize  'Nstat' -- can do this after loading the file

boxplot(Nstat[,1,], range=2)

## Log log scale
op <- sfsmisc::mult.fig(length(alphas))$old.par
for(i in seq_along(alphas)) {
    N.st <- Nstat[,i,]
    plot(NA,NA, xlim = range(V0s), ylim=range(N.st),
         main = substitute(N == no({"calls to_ " * rstable1()}) * " for  " * { alpha == A },
                           list(A = alphas[i])),
         xlab=expression(V[0]), ylab = "Nstat", type = "n", log="xy")
    boxplot(N.st, at = V0s, range=2, add = TRUE)
    lines(apply(N.st, 2, mean, trim=.05) ~ V0s, col = 2, lwd=2)
}
par(op)

str(Nst.mns <- apply(Nstat, 3:2, mean, trim=.05))
str(Nst.q95 <- apply(Nstat, 3:2, quantile, prob=.95))
str(Nst.max <- apply(Nstat, 3:2, max))
str(Nst.min <- apply(Nstat, 3:2, min))

## need  k * V0 :  what k
## in the mean:
rowMeans(Nst.mns / V0s) ## ?== Euler's  e == exp(1)
## maximally
rowMeans(Nst.max / V0s)

## is it  'e'  ??
plot(V0s, rowMeans(Nst.mns / V0s), log="x"); rug(V0s); abline(h=exp(1), col=2)
L <- V0s >= 20
plot(V0s[L], rowMeans(Nst.mns[L,] / V0s[L]), log="x"); rug(V0s[L]); abline(h=exp(1), col=2)
## well, as a limit probably

mtit <- substitute(N == `# `*group("{", rstable1() * " calls", "}") *
                   "   for  " * alpha %in% group("{", list(A,B,...,Z), "}"),
                   list(A=alphas[1], B=alphas[2], Z=tail(alphas,1)))

## normal scale
matplot (V0s, Nst.max, ylim=range(Nstat), col = "gray20",
         ylab="", main=mtit, type="l", lty=2); rug(V0s)
matlines(V0s, Nst.min, col = "gray20", lty=2)
abline(a=0,b=1, col = "tomato"); text(20,20, expression(y == x), col="tomato", adj=-.2)
matlines(V0s, Nst.mns, type = "b", lty=1)
matlines(V0s, Nst.q95, col = "gray80", lty=3)

## log log scale
matplot (V0s, Nst.max, ylim=range(Nstat), col = "gray20",
         ylab="", main=mtit, type="l", lty=2, log="xy"); rug(V0s)
matlines(V0s, Nst.min, col = "gray20", lty=2)
abline(a=0,b=1, col = "tomato"); text(20,20, expression(y == x), col="tomato", adj=-.2)
matlines(V0s, Nst.mns, type = "b", lty=1)
matlines(V0s, Nst.q95, col = "gray80", lty=3)


## or just the "residuals"
matplot(V0s, Nst.mns / V0s, type = "b", lty=1, log = "x")
rug(V0s)
