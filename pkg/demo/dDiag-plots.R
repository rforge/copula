p.dDiag <-
    function(family, d, n = 1000,
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
    x <- seq(0,1, length = n)
    yMat <- sapply(th., function(theta)
                   dDiag(x, cop=onacopulaL(family, list(theta, 1:d))))
    mainTit <-
        paste("Diagonal densities of ", family, "\n d = ",d, sep="")
    non.fin <- !is.finite(yMat)
    nF.mat <-
        if(any(has.nF <- apply(non.fin, 2, any))) {
            i.non.fin <- apply(non.fin, 2, which.max)
            i.non.fin[!has.nF] <- NA
            cat(mainTit,":\n non-finite values -- for (theta,  x >= *):\n")
            print(cbind(theta = th.[has.nF], x = x[i.non.fin[has.nF]]))
        }
    matplot(x, yMat, type = "l", col = cols, lty = 1,
            ylab="dDiag(x, *)", main= mainTit)
    if(any(has.nF)) { ## mark the end points
        i1 <- i.non.fin - 1
        points(x[i1], yMat[cbind(i1, seq_along(tau))],
               col = cols, cex = 1.5, lwd=2, pch = 4)
    }
    abline(h=0, lty=3)
    c.tau <- format(tau)
    lleg <- lapply(seq_along(tau), function(j) {
        cc <- if(has.nF[j]) {
            i <- i.non.fin[j]
            sprintf("  f(%.3g) = %g", x[i], yMat[i,j])
        } else ""
	substitute(list(tau == TAU, theta == TH) ~ COMM,
                   list(TAU= c.tau[j], TH= th.[j], COMM= cc))

    })
    legend("topleft", do.call(expression, lleg),
	   lty = 1, col=cols, bty="n")
    invisible(list(d = d, tau=tau, theta=th., x = x, y.dDiag = yMat,
                   has.non.fin = has.nF, non.fin.mat = nF.mat))
}

par(ask = dev.interactive(orNone = TRUE))

r.dDiag.3 <- lapply(nacopula:::c_longNames,
                    function(family) p.dDiag(family, d = 3))

r.dDiag.15 <- lapply(nacopula:::c_longNames,
                     function(family) p.dDiag(family, d = 15))

r.dDiag.75 <- lapply(nacopula:::c_longNames,
                     function(family) p.dDiag(family, d = 75))

