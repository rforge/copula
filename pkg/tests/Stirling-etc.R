library(nacopula)

### ---------------- Stirling numbers of the 1st kind ---------------------------

S1.10 <- c(0, -362880, 1026576, -1172700, 723680,
           -269325, 63273, -9450, 870, -45, 1)
stopifnot(sapply(0:10, Stirling1, n=10) == S1.10,
          Stirling1.all(10) == S1.10[-1])

options(str = strOptions(vec.len = 10, digits.d = 20)) # for ls.str() below

ls.str(nacopula:::.nacopEnv)
system.time(S  <- Stirling1(30, 7))# not zero
system.time(S. <- Stirling1(30, 7))# should be zero
stopifnot(identical(S, S.))

ls.str(nacopula:::.nacopEnv)

system.time(s1c <- Stirling1(100,10))
s1c
stopifnot(print(system.time(for(i in 1:20) S. <- Stirling1(100, 10))[[1]]) <= 0.010,
          identical(S., s1c))
system.time(s2c <- Stirling1(200,190)); s2c
stopifnot(print(system.time(for(i in 1:20) S. <- Stirling1(200,190))[[1]]) <= 0.010,
          identical(S., s2c))


### ---------------- Stirling numbers of the 2nd kind ---------------------------

S2.10 <- c(0, 1, 511, 9330, 34105, 42525, 22827, 5880, 750, 45, 1)
stopifnot(sapply(0:10, Stirling2, n=10, method="direct") == S2.10,
          sapply(0:10, Stirling2, n=10, method="lookup") == S2.10,
          Stirling2.all(10) == S2.10[-1])

ls.str(nacopula:::.nacopEnv)
system.time(S  <- Stirling2(30, 7))# not zero
system.time(S. <- Stirling2(30, 7))# should be zero
stopifnot(identical(S, S.),
          all.equal(S, Stirling2(30,7, method="direct"), tol=1e-15))

ls.str(nacopula:::.nacopEnv)

rbind(C.direct = system.time(Sd <- Stirling2(100,10, method="direct")),
      C.lookup = system.time(Sl <- Stirling2(100,10, method="lookup")))
## should be equal; and lookup time should be "zero" when called again:
stopifnot(all.equal(Sd, Sl, tol = 1e-15),
          print(system.time(for(i in 1:20) S. <- Stirling2(100, 10))[[1]]) <= 0.010)

## Here, the direct method already overflows, but the "lookup" still works
rbind(C.direct = system.time(Sd <- Stirling2(200,190, method="direct")),
      C.lookup = system.time(Sl <- Stirling2(200,190, method="lookup")))
Sd ; Sl
stopifnot(print(system.time(for(i in 1:20) S. <- Stirling2(200,190))[[1]]) <= 0.010)


### ---------------- Polylogarithm Function -------------------------------------

EQ <- function(x,y, tol = 1e-15) all.equal(x,y, tol=tol)

x <- (0:127)/128 # < 1
stopifnot(EQ(polylog(s =  1,  x, n.sum=10000), -log(1-x)),
	  EQ(polylog(s = -1, .1, n.sum=	 100), 10/81),
	  EQ(polylog(s = -1, .1, "neg"),       10/81),
	  EQ(polylog(x, -1, "neg"), x /(1-x)^2),
	  EQ(polylog(x, -2, "neg"), x*(1+x)/(1-x)^3),
	  EQ(polylog(x, -4, "neg"), x*(1+x)*(1+x*(10+x)) / (1-x)^5),
	  identical(	      polylog	  (x, -4, "neg"),
		    Vectorize(polylog,"z")(x, -4, "neg")),
	  identical(	      polylog	  (x, -4, "sum", n.sum=10000),
		    Vectorize(polylog,"z")(x, -4, "sum", n.sum=10000)),
          TRUE)


## Plots: ---
pdf("polylog-ex.pdf")

p.Li <- function(s.set, from = -2.6, to = 1/4, ylim = c(-1, 0.5),
                 colors = c("orange","brown", palette()), n = 201, ...)
{
    s.set <- sort(s.set, decreasing = TRUE)
    s <- s.set[1] # <_ for auto-ylab
    curve(polylog(x, s, method="neg"), from, to,
          col=colors[1], ylim=ylim, n=n, ...)
    abline(h=0,v=0, col="gray")
    for(is in seq_along(s.set)[-1])
        curve(polylog(x, s=s.set[is], method="neg"), add=TRUE, col = colors[is], n=n)
    s <- rev(s.set)
    legend("bottomright", paste("s =", s), col=colors[2-s], lty=1, bty="n")
}

## yellow is unbearable (on white):
palette(local({p <- palette(); p[p=="yellow"] <- "goldenrod"; p}))

## Wikipedia page plot (+/-):
p.Li(1:-3, ylim= c(-.8, 0.6), colors = c(2:4,6:7))

## and a bit more:
p.Li(1:-5)

## For the range we need it:
ccol <- c(NA,NA, rep(palette(),10))
p.Li(-1:-20, from=0, to=.99, colors=ccol, ylim = c(0, 10))
## log-y scale:
p.Li(-1:-20, from=0, to=.99, colors=ccol, ylim = c(.01, 1e7), log="y")
