#include <Rmath.h>

#include "nacopula.h"

/* From: Marius Hofert <marius.hofert@math.ethz.ch>
 * Date: Mon, 26 Apr 2010 12:15:51 +0200 */

/* sample a Log(p) distribution, see Kemp (1981), algorithm ``LK'' : */
double rLog(double p) {

    /* NOTA BENE: the *caller* of this function *MUST* use
     * --------- GetRNGstate(); .... ; PutRNGstate();
     */

    if(p <= 0. ||  p >= 1.) {
	error("rLog(): p must be inside (0,1)");
	return -1.;/* -Wall */
    }
    else {
	double U=unif_rand();
	if(U > p) {
	    return 1.;
	}
	else {
	    double Q = - expm1(log1p(- p) * unif_rand());
		/* == 1. - exp(log1p(- p) * unif_rand())
		 * == 1. - pow(1. - p, unif_rand());
		 */
	    return(U < Q*Q
		   ? floor(1. + log(U)/log(Q))
		   : ((U > Q) ? 1. : 2.));
	}
    }
}

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
