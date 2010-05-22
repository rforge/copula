#include <Rmath.h>

#include "nacopula.h"

/**
 * Sample a Log(p) distribution with the algorithm "LK" of Kemp (1981).
 * Note: The caller of this function must use GetRNGstate() and PutRNGstate().
 * @param p in (0,1)
 * @return a random variate from Log(p)
 * @author Marius Hofert, Martin Maechler
*/
double rLog(double p) {
    if(p <= 0. ||  p >= 1.) {
	error("rLog(): p must be inside (0,1)");
	return -1.; /**< -Wall */
    }
    else {
	double U=unif_rand();
	if(U > p) {
	    return 1.;
	}
	else {
	    double Q = - expm1(log1p(- p) * unif_rand());
		/**
		 * == 1. - exp(log1p(- p) * unif_rand())
		 * == 1. - pow(1. - p, unif_rand())
		 */
	    return(U < Q*Q
		   ? floor(1. + log(U)/log(Q))
		   : ((U > Q) ? 1. : 2.));
	}
    }
}

/**
 * Generate a vector of variates from a Log(p) distribution.
 * @param n_ sample size
 * @param p_ parameter p in (0,1)
 * @return vector of random variates from Log(p)
 * @author Martin Maechler
*/
SEXP rLog_c(SEXP n_, SEXP p_) {
    int n = asInteger(n_);
    double p = asReal(p_);
    SEXP res = PROTECT(allocVector(REALSXP, n));
    double* X = REAL(res);

    GetRNGstate();

    for(int i=0; i < n; i++)
	X[i] = rLog(p);

    PutRNGstate();
    UNPROTECT(1);
    return res;
}
