#ifndef NACOPULA_DEFS_H
#define NACOPULA_DEFS_H

#include <R.h>
#include <Rinternals.h>

SEXP sinc(SEXP x_);
SEXP A_Zolotarev(SEXP x_, SEXP alpha, SEXP I_alpha);

SEXP rLog_c     (SEXP n, SEXP p);
SEXP rFJoe_c    (SEXP n, SEXP alpha);
SEXP rFFrank_c  (SEXP n, SEXP theta_0, SEXP theta_1);
SEXP retstable_c(SEXP V0, SEXP alpha);
SEXP rstable_b1 (SEXP n, SEXP alpha, SEXP gamma);

/* -- C API --- for "us" but maybe other R packages ---> need to "export" it via ../inst/include/ */
double R_sinc(double x);

double rFJoe(double alpha,
	     double iAlpha,     /* := 1 - alpha */
	     double gamma_1_a); /* == Gamma(1 - alpha) == Gamma(iALpha) */
void rFJoe_vec(double V[], const int n,
	       const double alpha, const double iAlpha /* := 1 - alpha */);

/* two versions for now : */
double S (const double alpha, const double gamma);
double S_(const double alpha, const double gamma);

void C_retstable  (double *V, const double V0[], double alpha, int n);
void C_retstable_R(double *V, const double V0[], double alpha, int n);

void rstable_b1_vec(double V[], const int n,
		    const double alpha, const double gamma);

double rLog(double p);

double rFFrank(double p, double alpha, double iAlpha /* == 1 - alpha */,
	       int theta_0_le_1);
void rFFrank_vec(double *X, const int n,
		 const double theta_0, const double theta_1);

#endif
