## Utility for  polyG()  in ./R/aux-acopula.R :
sign.fafac <- function(alpha, d) { ## also  on (d ; k:= 1:d ; k1:= 0:(d-1))
    stopifnot(length(alpha)==1, 0 < alpha, alpha < 1,
              length(d) == 1, d == round(d), d >= 1)
    k <- 1:d
    k1 <- k - 1L
    s <- unlist(lapply(alpha*k, function(z) prod(z-k1)))
    ss <- sign(s)
    sn <- (-1)^d * (2*(floor(alpha*k) %% 2) - 1)
    ##    --------------------------------------
    stopifnot(ss[s != 0] == sn[s != 0])
    ss
}

sign.fafac(0.8, 4)
sign.fafac(0.8, 5)
sign.fafac(0.8, 6)
sign.fafac(0.8, 16)
sign.fafac(0.8, 17)

## Ok: It is "proven"
##
##    sign(s) == (-1)^d * (2*(floor(alpha*k) %% 2) - 1)
##
## at whenever sign(s) != 0  {and for the s == 0  case, sign does not matter}
