#include <Rmath.h>

#include "nacopula.h"

/** rejection for F for Frank's family
 * @param p
 * @param alpha
 * @param theta0le1

 * @return
 */
double rFFrank(double p, double alpha, double iAlpha /* == 1 - alpha */, int theta_0_le_1)
{
    double U, V;
    if(theta_0_le_1) {
	do {
	    U = unif_rand();
	    V = rLog(p);
	} while (U*(V-alpha) > 1./beta(V, iAlpha));

    } else {
	double gamma_1_a = gammafn(iAlpha);
	do {
	    U = unif_rand();
	    V = rFJoe(alpha, iAlpha, gamma_1_a);
	} while(U > pow(p, V-1.));
    }
    return V;
}

void rFFrank_vec(double *V, const int n,
		 const double theta_0, const double theta_1)
{
    double p  =  - expm1(-theta_1),
	alpha = theta_0 / theta_1, iAlpha = (theta_1 - theta_0) / theta_1;
    int th_0_le_1 = (theta_0 <= theta_1);

    if(n >= 1) {
	GetRNGstate();

	for(int i=0; i < n; i++)
	    V[i] = rFFrank(p, alpha, iAlpha, th_0_le_1);

	PutRNGstate();
    }
    return;
}

SEXP rFFrank_c(SEXP n_, SEXP theta_0_, SEXP theta_1_) {
    int n = asInteger(n_);
    double theta_0 = asReal(theta_0_),
	   theta_1 = asReal(theta_1_);
    SEXP res = PROTECT(allocVector(REALSXP, n));
    if(n >= 1)
	rFFrank_vec(REAL(res), n, theta_0, theta_1);
    UNPROTECT(1);
    return res;
}
