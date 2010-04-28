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

/* MM: this is really not (yet) useful:- CPU time saving is close to 0
 * -- if you had used  Rprof() as I said,  you'd have seen that really only
 *  rstable1() is the costly part ! */

/**
* <description>
* @title Sample  length(V0)  variates from exponentially-tilted Stable distribution
*
* <details> Sample from
* S~(alpha,1, (cos(alpha*pi/2)V0/m)^(1/alpha), V0*I(alpha==1), I(alpha!=1);1)
* with corresponding Laplace-Stieltjes transform
* exp(-V0*((1+t)^alpha-1))

* @param V0_  numeric(n) {R vector}
* @param alpha parameter {R numeric(1)}
* @return R numeric(n) vector
*/
SEXP retstable_c(SEXP V0_, SEXP alpha)
{
    int n = LENGTH(PROTECT(V0_ = coerceVector(V0_, REALSXP)));
    double alp = asReal(alpha);
    SEXP V_ = allocVector(REALSXP, n); /* the result */
    PROTECT(V_);

    C_retstable(REAL(V_), REAL(V0_), alp, n);
    /* or
       C_retstable_R(REAL(V_), REAL(V0_), alp, n);
    */

    UNPROTECT(2);
    return(V_);
}


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

/* computes a stable  S(alpha, beta=1, gamma, \I_{\{\alpha=1\}}; 1) distributed random variate
 * with exponent alpha, skewness parameter beta=1, scale parameter gamma, and
 * shift parameter delta = I_{\alpha = 1} --- adapted from gsl_ran_levy_skew */
double S(const double alpha, const double gamma)
{
    double delta=0.;
    double Theta=M_PI*(unif_rand() - 0.5), W;
    do {
	/* computes an variate from an Exp(lambda = 1) - distribution: */
	W = -log(unif_rand());
    } while (W == 0.); /* U = unif_rand() == 1  should never happen */

    if(alpha == 1.) { /* NOTA BENE: Marius had this wrong .. */
#define p2 M_PI_2
	double p2bT = p2 + /* beta * */ Theta;
	return (p2bT*tan(Theta) - /* beta * */ log((p2 * W * cos(Theta))/p2bT))/p2;
#undef p2
    }
    else{
	double t = /* beta * */ tan(M_PI_2*alpha);
	double th0 = atan(t)/alpha;
	double s = pow(1.+t*t, 1./(2.*alpha));
	double x = s * sin(alpha*(Theta+th0)) / pow(cos(Theta), 1./alpha) *
	    pow(cos(Theta-alpha*(Theta+th0)) / W, (1.-alpha)/alpha);
	return delta + gamma*x;
    }
}


void C_retstable(double *V, const double V0[], double alpha, int n)
{
    double c0 = cos(M_PI_2 * alpha);

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
	double gamma = pow(c0*V0[i]/m, 1./alpha);
	for(int k= 0; k < m; k++) {
	    double V01_k, U;
	    /* standard rejection */
	    do {
		/* sample from the distribution corresponding to the Laplace-Stieltjes
		   transform exp(-(V0_0/m)*t^alpha) */
		V01_k = S(alpha, gamma); /* rstable(1, alpha, beta=1, gamma) : */
		U = unif_rand();

	    } while (U > exp(- V01_k));

	    V[i] += V01_k; /* V__ := sum_{k'=1}^k  X  */
	}
    } /* end for(i .. n) */
    PutRNGstate();
    return;
}
