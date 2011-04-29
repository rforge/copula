require(nacopula)

### We have other implicit estimation "tests" in  demo(estimation.gof)
## i.e., ../demo/estimation.gof.R
##       ~~~~~~~~~~~~~~~~~~~~~~~~

## Frank & mle ---- log-density is *REALLY* not good enough:

(cop <- onacopulaL("Frank",list(NA, 1:5)))
uu <- rbind(rep(.987, 5))
th <- seq(1,60, by=1/4)
## log-likelihood:
l.th <- sapply(th, cop@copula@dacopula, u = uu, log=TRUE,
## this was default ~ 2011-03:
               method = "negI-s-Stirling", Li.log.arg=FALSE)
plot(l.th ~ th, type = "l",   ## jumps from a noisy "12" to +Inf --- the horror
     col=2, xlab = expression(theta))
l.th2 <- sapply(th, cop@copula@dacopula, u = uu, log=TRUE,
               method = "negI-s-Stirling", Li.log.arg=TRUE)
lines(l.th2 ~ th, col="green2")

## using  MC  (instead of polylog(.))  also ``works'' {the function is noisy ..}:
l.th.MC <- sapply(th, cop@copula@dacopula, u = uu, n.MC = 1e4, log=TRUE)
matplot(th, cbind(l.th, l.th2, l.th.MC), type = "l",
        col=c("red", "green2", "blue"), lty=c(2,1,3))

cbind(th, l.th, l.th2)
