library(nacopula)


##==== Original "Tests" from Marius ====

set.seed(1) # V0() use random numbers

copAMH@psi(1, 0.5)
copAMH@psiInv(0.5, 0.5)
copAMH@nestConstr(0.5, 0.2)
copAMH@V0(10, 0.5)
copAMH@V01(10, 0.5, 0.8)
copAMH@tau(0.5)
copAMH@tauInv(0.2)
copAMH@lTDC(0.5)

copClayton@psi(1, 0.5)
copClayton@psiInv(0.5, 0.5)
copClayton@nestConstr(0.2, 0.5)
(V0 <- copClayton@V0(10, 0.5))
copClayton@V01(V0,  0.5,  0.8)
copClayton@tau(0.5)
copClayton@tauInv(0.2)
copClayton@lTDC(0.5)
copClayton@lTDCInv(0.5)

n <- 1000
(theta0 <- copClayton@tauInv(0.05))
(theta1 <- copClayton@tauInv(0.1))
(V0 <- copClayton@V0(n, theta0)) ## fast
startclock <- proc.time()
cat('Time elapsed: ', startclock,'\n')
V01 <- copClayton@V01(V0, theta0, theta1) ## slow
runtime <- proc.time()-startclock
V01
cat('Time elapsed: ', runtime,'\n') # for ``statistical reasons''
