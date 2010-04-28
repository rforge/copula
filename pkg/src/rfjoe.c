#include <Rmath.h>

#include "nacopula.h"

double rFJoe(double alpha,
	     double iAlpha,    /* := 1 - alpha */
	     double gamma_1_a) /* == Gamma(1 - alpha) == Gamma(iALpha) */
{
    double U, I_al = 1./alpha;

    /* NOTA BENE: the *caller* of this function *MUST* use
     * --------- GetRNGstate(); .... ; PutRNGstate();
     */

    /* FIXME(MM): (for alpha not too close to 1): re-express using 1-U !*/
    U = unif_rand();
    if(U <= alpha)
	return 1.;
    else { /* alpha < U < 1 */
	double Ginv = pow((1-U)*gamma_1_a, -I_al);
	double fGinv = floor(Ginv);
	if(1-U < 1./(fGinv*beta(fGinv, iAlpha)))
	    return ceil(Ginv);
	else return fGinv;
    }
}

void rFJoe_vec(double V[], const int n,
	       const double alpha, const double iAlpha /* := 1 - alpha */)
{
    if(n >= 1) {
	double G1_a = gammafn(iAlpha);
	GetRNGstate();

	for(int i=0; i < n; i++)
	    V[i] = rFJoe(alpha, iAlpha, G1_a);

	PutRNGstate();
    }
    return;
}


/**
* <description>
*
* <details>
* @title Sample F with probability mass function p_k=\binom{alpha}{k}(-1)^{k-1},
*   k in IN, for Joe's family
* 	Note: should be *fast* as it is used as building block in many places
* @param n sample size
* @param alpha parameter
* @return numeric(n) vector
*/
SEXP rFJoe_c(SEXP n, SEXP alpha)
{
    int nn = asInteger(n);
    double alp = asReal(alpha);
    SEXP res = PROTECT(allocVector(REALSXP, nn));

    rFJoe_vec(REAL(res), nn, alp, 1. - alp);

    UNPROTECT(1);
    return(res);
}
