library(nacopula)

myCop <- copAMH
theta(myCop) <- 0.5

cop2 <- setTheta(copAMH, value = 0.5) # is maybe more natural
stopifnot( identical(cop2, myCop) )

setGeneric("psi", function(cop) standardGeneric("psi"))
setMethod(psi, "ACopula",
          function(cop) {function(t) cop@psi(t, theta = cop@theta)})
psi(myCop) ## is a *function*
psi(myCop)(0:4)
curve(psi(myCop)(x), 0, 4)
## but this can also be done directly [ => same curve "on top" :]
curve(myCop@psi(x, theta = myCop@theta),  0, 4, col = 2, add=TRUE)

## Kendall's tau --
(ii <- copClayton@paraInterval) # upper bound Inf
th <- seq(ii[1], 50, length = 201)
plot(th, copClayton@tau(th), type = "l",
     main = expression(tau["Clayton"](theta)))

p.Tau <- function(cop, n = 201, xlim = pmin(paraI, 50), ...) {
    stopifnot(is(cop, "ACopula"))
    paraI <- cop@paraInterval
    theta <- seq(xlim[1], xlim[2], length.out = n)
    tit <- substitute(tau[NAME](theta), list(NAME = cop@name))
    plot(theta, cop@tau(theta), type = "l", main = tit, ...)
}

p.Tau(copAMH)
p.Tau(copClayton)
if(FALSE) ## FIXME
p.Tau(copFrank)# -> error in integrate !
p.Tau(copGumbel)
if(FALSE) ## FIXME
p.Tau(copJoe)# another error

