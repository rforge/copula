library(nacopula)

S1.10 <- c(0, -362880, 1026576, -1172700, 723680,
           -269325, 63273, -9450, 870, -45, 1)
stopifnot(sapply(0:10, Stirling1, n=10) == S1.10)

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


###  --------------------- Stirling  2nd kind ------------------------------

S2.10 <- c(0, 1, 511, 9330, 34105, 42525, 22827, 5880, 750, 45, 1)
stopifnot(sapply(0:10, Stirling2, n=10, method="direct") == S2.10,
          sapply(0:10, Stirling2, n=10, method="lookup") == S2.10)

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


