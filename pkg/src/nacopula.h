#ifndef NACOPULA_DEFS_H
#define NACOPULA_DEFS_H

#include <R.h>
#include <Rinternals.h>

SEXP sinc_c(SEXP x_);
SEXP A__c(SEXP x_, SEXP alpha, SEXP I_alpha);

SEXP rstable_c(SEXP n, SEXP alpha, SEXP gamma);
SEXP retstable_MH_c(SEXP V0_, SEXP h, SEXP alpha);
SEXP retstable_LD_c(SEXP V0_, SEXP h, SEXP alpha);

SEXP rLog_c     (SEXP n, SEXP p);
SEXP rFJoe_c    (SEXP n, SEXP alpha);
SEXP rFFrank_c  (SEXP n, SEXP theta_0, SEXP theta_1);

/* -- C API --- for "us" but maybe other R packages ---> need to "export" it via ../inst/include/ */
double sinc_MM(double x);
double A_(double x, double alpha);
double BdB0(double x, double alpha);

double rstable0(double alpha);
double rstable (double alpha, double gamma);
void rstable_vec(double S[], const int n,
		 const double alpha, const double gamma);
void retstable_MH(double *St, const double V0[], double h, double alpha, int n);
void retstable_LD(double *St, const double V0[], double h, double alpha, int n);

double rLog(double p);
double rFJoe(double alpha,
	     double iAlpha,     /* := 1 - alpha */
	     double gamma_1_a); /* == Gamma(1 - alpha) == Gamma(iALpha) */
void rFJoe_vec(double V[], const int n,
	       const double alpha, const double iAlpha /* := 1 - alpha */);
double rFFrank(double p, double alpha, double iAlpha /* == 1 - alpha */,
	       int theta_0_le_1);
void rFFrank_vec(double *X, const int n,
		 const double theta_0, const double theta_1);

#endif
