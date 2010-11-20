/*
 Copyright (C) 2010 Marius Hofert and Martin Maechler

 This program is free software; you can redistribute it and/or modify it under
 the terms of the GNU General Public License as published by the Free Software
 Foundation; either version 3 of the License, or (at your option) any later
 version.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 details.

 You should have received a copy of the GNU General Public License along with
 this program; if not, see <http://www.gnu.org/licenses/>.
*/

#ifndef NACOPULA_DEFS_H
#define NACOPULA_DEFS_H

#include <R.h>
#include <Rinternals.h>

/**< For internationalized messages */
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("Matrix", String)
#else
#define _(String) (String)
#define dngettext(pkg, String, StringP, N) (N > 1 ? StringP : String)
#endif

SEXP sinc_c(SEXP x_);
SEXP A__c(SEXP x_, SEXP alpha, SEXP I_alpha);

SEXP rstable_c(SEXP n, SEXP alpha);
SEXP retstable_c(SEXP V0_, SEXP h, SEXP alpha, SEXP method);

SEXP rLog_c     (SEXP n, SEXP p, SEXP Ip);
SEXP rFJoe_c    (SEXP n, SEXP alpha);
SEXP rFFrank_c  (SEXP n, SEXP theta_0, SEXP theta_1);

/**
 * C API---for "us" but maybe also other R packages
 * "export" it via ../inst/include/
*/
double sinc_MM(double x);
double A_(double x, double alpha);
double BdB0(double x, double alpha);

double rstable0(double alpha);
double rstable (double alpha);
void rstable_vec(double S[], const int n, const double alpha);
void retstable_MH(double *St, const double V0[], double h, double alpha, int n);
void retstable_LD(double *St, const double V0[], double h, double alpha, int n);

double rLog(double p, double Ip);
double rFJoe(double alpha,
	     double iAlpha,     /* := 1 - alpha */
	     double gamma_1_a); /* == Gamma(1 - alpha) == Gamma(iALpha) */
void rFJoe_vec(double V[], const int n,
	       const double alpha, const double iAlpha /* := 1 - alpha */);
double rFFrank(double p, double Ip, /* = 1 - p */
	       double alpha, double iAlpha /* := 1 - alpha */,
	       int theta_0_le_1);
void rFFrank_vec(double *X, const int n,
		 const double theta_0, const double theta_1);

#endif
