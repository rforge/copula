p.dDiag <-
    function(family, d, n = 1000, log = FALSE,
             tau = c(0.01,
                     ## First tau = 0.01 for all families
                     if(family=="AMH") c((1:3)/10, 0.33) else
                     c((1:9)/10, .95, 0.99, .999)),
             cols = adjustcolor(colorRampPalette(c("red", "orange", "blue"),
                                space = "Lab")(length(tau)), 0.8))
{
    stopifnot(length(d) == 1, d == as.integer(d), d >= 2, n >= 10)
    cop <- getAcop(family)
    th. <- cop@tauInv(tau)
    u <- seq(0,1, length = n)
    mainTit <-
        paste("Diagonal densities of ", family, "\n d = ",d,
              if(log)", log = TRUE", sep="")
    c.tau <- format(tau)
    c.th  <- format(th., scientific=FALSE, digits = 4)

    ## yMat <- sapply(th., function(theta)
    ##                dDiag(u, cop=onacopulaL(family, list(theta, 1:d)),
    ##                      log=log))
    ## The "generic" Archimedean  dDiag():
    yM.a <- sapply(th., function(theta)
                   nacopula:::dDiagA(u, d=d, cop= setTheta(cop, theta), log=log))
    ## The dDiag-*s*lot :
    yM.s <- sapply(th., function(theta)
                   cop@dDiag(u, theta=theta, d=d, log=log))

    p.oneMat <- function(yMat, label) {
        thisTit <- paste(label, mainTit, sep=": ")
        non.fin <- !is.finite(yMat)
        ## Rather: in log-scale, -Inf  may be ok :
        if(log) non.fin <- non.fin & (is.na(yMat) | yMat != -Inf)
        nF.mat <-
            if(any(has.nF <- apply(non.fin, 2, any))) {
                i.non.fin <- apply(non.fin, 2, which.max)
                i.non.fin[!has.nF] <- NA
                cat(thisTit,":\n non-finite values -- for (theta,  u >= *):\n")
                print(cbind(theta = th.[has.nF], u = u[i.non.fin[has.nF]]))
            }
        matplot(u, yMat, type = "l", col = cols, lty = 1,
                ylab="dDiag(u, *)", main= thisTit)
        if(any(has.nF)) { ## mark the end points
            i1 <- i.non.fin - 1
            points(u[i1], yMat[cbind(i1, seq_along(tau))],
                   col = cols, cex = 1.5, lwd=2, pch = 4)
        }
        abline(h=0, lty=3)
        lleg <- lapply(seq_along(tau), function(j) {
            cc <- if(has.nF[j]) {
                i <- i.non.fin[j]
                sprintf("  f(%.3g) = %g", u[i], yMat[i,j])
            }
            substitute(list(tau == TAU, theta == TH) ~ COMM,
                       list(TAU= c.tau[j], TH= c.th[j], COMM= cc))

        })
        legend(if(log)"bottomright" else "topleft",
               do.call(expression, lleg),
               lty = 1, col=cols, bty="n")
    }## end{ p.oneMat() }

    p.oneMat(yM.a, "dDiagA()")
    p.oneMat(yM.s, "cop @ dDiag()")

    invisible(list(d = d, tau=tau, theta=th., u = u, dDiag.a = yM.a, dDiag.s = yM.s))
}

par(ask = dev.interactive(orNone = TRUE))

r.dDiag.3 <- lapply(nacopula:::c_longNames,
                    function(family) p.dDiag(family, d = 3))
r.dDiag.4.L <- lapply(nacopula:::c_longNames,
                    function(family) p.dDiag(family, d = 4, log=TRUE))

r.dDiag.15 <- lapply(nacopula:::c_longNames,
                     function(family) p.dDiag(family, d = 15))

r.dDiag.75 <- lapply(nacopula:::c_longNames,
                     function(family) p.dDiag(family, d = 75))

r.dDiag.200.L <- lapply(nacopula:::c_longNames,
                    function(family) p.dDiag(family, d = 200, log=TRUE))
