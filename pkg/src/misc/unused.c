#include <Rmath.h>

#include "nacopula.h"

/* lang5 will probably be defined in future versions of R :*/
#ifndef lang5
static SEXP lang5(SEXP s, SEXP t, SEXP u, SEXP v, SEXP w)
{
    PROTECT(s);
    s = LCONS(s, list4(t, u, v, w));
    UNPROTECT(1);
    return s;
}
#endif


/* MM: this is really not useful, as the time is spent in the call to R's
 *     rstable1() !
 * keep it as example / for didactical reasons */
void C_retstable_R(double *V, const double V0[], double alpha, int n)
{
    double c0 = cos(M_PI_2 * alpha);
    SEXP one, alpha_R;
    one     = PROTECT(ScalarReal(1.));
    alpha_R = PROTECT(ScalarReal(alpha));

    GetRNGstate();
    for(int i = 0; i < n; i++) {
	/* find m := optimal number of summands -- now using asymptotic formula :*/
	double m = fmax(1., floor(0.5 + V0[i]));
	/* Previously:
	 * double m;
	 * if(V0[i] <= 1)
	 *     m= 1.;
	 * else {
	 *     double fV= floor(V0[i]), cV= ceil(V0[i]);
	 *     if(pow(exp(-V0[i]), 1./cV - 1./fV) <= cV/fV)
	 * 	m= fV;
	 *     else
	 * 	m= cV;
	 * } */

	/* apply standard rejection m times */
	SEXP gamma_, rstableCall;
	PROTECT(gamma_ = ScalarReal(pow(c0*V0[i]/m, 1./alpha)));
	PROTECT(rstableCall = lang5(install("rstable1"), one, /* alpha=*/ alpha_R,
				    /* beta=*/ one, /* gamma= */ gamma_));
	for(int k= 0; k < m; k++) {
	    double V01_k, U;
	    /* standard rejection */
	    do {
		/* sample from the distribution corresponding to the Laplace-Stieltjes
		   transform exp(-(V0_0/m)*t^alpha) */
		V01_k = asReal(eval(rstableCall, R_GlobalEnv));
		U = unif_rand();

	    } while (U > exp(- V01_k));

	    V[i] += V01_k; /* V__ := sum_{k'=1}^k  X  */
	}
	UNPROTECT(2);
    } /* end for(i .. n) */
    PutRNGstate();
    UNPROTECT(1);
    return;
}

