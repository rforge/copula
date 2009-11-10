library(nacopula)


##==== Original "Tests" from Marius ====

set.seed(1) # V0() use random numbers

copAMH@psi(1, 0.5)
psiinvA(0.5, 0.5)
nconstrA(0.5, 0.2)
V0A(10, 0.5)
V01A(10, 0.5, 0.8)
tauA(0.5)
tauinvA(0.2)
lambdalA(0.5)

psiC(1, 0.5)
psiinvC(0.5, 0.5)
nconstrC(0.2, 0.5)
(V0 <- V0C(10, 0.5))
V01C(V0,  0.5,  0.8)
tauC(0.5)
tauinvC(0.2)
lambdalC(0.5)
lambdalinvC(0.5)

n <- 1000
(theta0 <- tauinvC(0.05))
(theta1 <- tauinvC(0.1))
(V0 <- V0C(n, theta0)) ## fast
startclock <- proc.time()
cat('Time elapsed: ', startclock,'\n')
V01 <- V01C(V0, theta0, theta1) ## slow
runtime <- proc.time()-startclock
V01
cat('Time elapsed: ', runtime,'\n') # for ``statistical reasons''
