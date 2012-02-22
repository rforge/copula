#ifndef COPULA_DEFS_H
#define COPULA_DEFS_H

#include <R.h>

#include "Anfun.h"
#include "gof.h"
#include "set.utils.h"
#include "empcop.stat.h"

// ./logseries.c : __FIXME__ also have rLog_vec_c()  from nacopula
void rlogseries_R(int *n, double *alpha, int *val);

// ./fgm.c:
void validity_fgm(int *p, double *alpha, int *valid);
void rfgm(int *p, double *alpha, int *n, double *x);

// ./evtest.c : ------------------------------------------------
void evtest(double *U, int *n, int *p, double *g, int *m,
	    int *N, double *tg, int *nt, double *s0, int *der2n,
	    double *o, double *stat);

void evtestA(double *U, double *V, int *n, double *u, double *v,
	     int *m, int *CFG, int *N, double *s0);

void evtestA_derA(double *U, double *V, int *n, double *u, double *v,
		  int *m, int *CFG, int *N, double *s0);

void evtestA_stat(double *U, double *V, int *n, double *u, double *v, int *m,
		  int *CFG, double *stat, double *offset);

// "_C": nameclash ...
void evTestAA_C(double *U, double *V, int *n, double *t, int *m,
		int *N, double *s0);
void evTestAA_derA(double *U, double *V, int *n, double *t, int *m,
		   int *N, double *s0);
void evTestAA_stat(double *S, double *T, int *n, double *t, int *m,
		   double *stat);

// ./exchtest.c : ------------------------------------------------
void evsymtest(double *U, double *V, int *n, double *t, int *m,
	       int *CFG, int *N, double *s0);

void evsymtest_derA(double *U, double *V, int *n, double *t, int *m,
		    int *CFG, int *N, double *s0);

void evsymtest_stat(double *S, double *T, int *n, double *t, int *m,
		    int *CFG, double *stat);

void exchtestCn(double *U, double *V, int *n, double *u, double *v,
		int *m, int *N, double *s0);

void exchtestCn_stat(double *U, double *V, int *n, double *u, double *v,
		     int *m, double *stat);

// R_debye.c : -----------------------------------------------------------------
// "_C": nameclash - already have R level 'debye_1'
void debye_1_C(double *x, int *len, double *val, double *err, int *status);
void debye_2(double *x, int *len, double *val, double *err, int *status);
void debye_3(double *x, int *len, double *val, double *err, int *status);
void debye_4(double *x, int *len, double *val, double *err, int *status);


#endif
