/*
  Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan

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

/**
 * @file   gof.c
 * @author Ivan Kojadinovic
 * @date   2012
 *
 * @brief  Some goodness-of-fit tests for copulas and EV copulas
 *         based on the parametric bootstrap and on a multiplier
 *         bootstrap
 *
 */


#include <R.h>
#include <Rmath.h>

#include "Anfun.h"
#include "gof.h"
#include "empcop.h"

/**
 * Cramer-von Mises test statistic
 *
 * @param n sample size
 * @param p dimension
 * @param U pseudo-observations
 * @param Ctheta values of the fitted copula at U
 * @param stat value of the test statistic
 * @author Ivan Kojadinovic
 */
void cramer_vonMises(int *n, int *p, double *U, double *Ctheta,
		     double *stat)
{
  double s = 0.;
  for(int k=0; k < *n; k++) {
    double diff = multCn(U,*n,*p,U,*n,k,0.0) - Ctheta[k];
    s += diff * diff;
  }
  *stat = s;
}

/**
 * Cramer-von Mises test statistic (grid version)
 * Used mainly for testing -- Not really needed anymore
 *
 * @param p dimension
 * @param U pseudo-observations
 * @param n sample size
 * @param V grid
 * @param m grid size
 * @param Ctheta values of the fitted copula at V
 * @param stat value of the test statistic
 * @author Ivan Kojadinovic
 */
void cramer_vonMises_grid(int *p, double *U, int *n, double *V, int *m,
			  double *Ctheta, double *stat)
{
  double s = 0.;
  for (int k=0; k < *m; k++) {
    double diff = multCn(U,*n,*p,V,*m,k,0.0) - Ctheta[k];
    s +=  diff * diff;
  }
  *stat = s * (*n) / (*m);
}

/**
 * Multiplier bootstrap for the GoF testing
 * TODO: rewrite so that n is the sample size!
 *
 * @param p dimension
 * @param u0 pseudo-observations (m x p)
 * @param m sample size
 * @param u grid (n x p)
 * @param n grid size
 * @param influ influence matrix (n x m)
 * @param N number of multiplie replications
 * @param s0 N replications of the test statistic
 * @author Ivan Kojadinovic
 */
void multiplier(int *p, double *u0, int *m, double *u, int *n,
		double *influ, int *N, double *s0)
{
  int i, j, k, l, ind;
  double *influ_mat = Calloc((*m) * (*n), double);
  double *random = Calloc(*m, double);
  double *v1 = Calloc(*p, double);
  double *v2 = Calloc(*p, double);
  double *der = Calloc(*p, double);
  double mean, process, invsqrtm = 1.0/sqrt(*m);

  /* influence matrix */
  for (j = 0; j < *n; j++) /* loop over the grid points */
    {
      /* derivatives wrt args */
      for (k = 0; k < *p; k++)
	{
	  v1[k] = u[j + k * (*n)];
	  v2[k] = v1[k];
	}
      for (k = 0; k < *p; k++)
	{
	  v1[k] += invsqrtm;
	  v2[k] -= invsqrtm;
	  der[k] = der_multCn(u0, *m, *p, v1, v2, 2 * invsqrtm);
	  v1[k] -= invsqrtm;
	  v2[k] += invsqrtm;
	}

      for (i = 0; i < *m; i++) /* loop over the data */
	{
	  influ_mat[i + j * (*m)] = 0.0;
	  ind = 1;
	  for (k = 0; k < *p; k++)
	    {
	      ind *= (u0[i + k * (*m)] <= u[j + k * (*n)]);
	      influ_mat[i + j * (*m)] -= der[k] * (u0[i + k * (*m)] <= u[j + k * (*n)]);
	    }
	  influ_mat[i + j * (*m)] += ind; /* - influ[j + i * (*n)];*/
	  influ[j + i * (*n)] *= invsqrtm;
	  influ_mat[i + j * (*m)] *= invsqrtm;
	}
    }

  GetRNGstate();

  /* generate N approximate realizations */
  for (l=0;l<*N;l++)
    {
      /* generate m variates */
      mean = 0.0;
      for (i=0;i<*m;i++)
	{
	  random[i] = norm_rand(); /*(unif_rand() < 0.5) ? -1.0 : 1.0 ;*/
	  mean += random[i];
	}
      mean /= *m;

      /* realization number l */
      s0[l] = 0.0;
      for (j=0;j<*n;j++)
	{
	  process = 0.0;
	  for (i=0;i<*m;i++)
	    process += (random[i] - mean) * influ_mat[i + j * (*m)]
	      - random[i] * influ[j + i * (*n)];
	  s0[l] += process * process;
	}
      s0[l] /= *n;
    }

  PutRNGstate();

  Free(influ_mat);
  Free(random);
  Free(v1);
  Free(v2);
  Free(der);
}

/// Goodness-of-fit tests for extreme-value copulas

/**
 * Cramer-von Mises test statistic based on the Pickands estimator
 * used in the GOF for EV copulas
 * stat = n/m \sum_{i=0}^{m-1} [An(i/m) - A_{theta_n}(i/m)]^2
 *
 * @param n sample size
 * @param m grid size
 * @param S unit Fréchet pseudo-observations
 * @param T unit Fréchet pseudo-observations
 * @param Atheta values of the fitted A at the grid points
 * @param stat value of the test statistic
 * @author Ivan Kojadinovic
 */
void cramer_vonMises_Pickands(int *n, int *m, double *S,
			      double *T, double *Atheta,
			      double *stat)
{
  int i;
  double t, Ac, Au, dc, du,
    invA0 = inv_A_Pickands(*n, S, T, 0.0),
    invA1 = inv_A_Pickands(*n, S, T, 1.0);

  stat[0] = 0.0; stat[1] = 0.0;
  for (i=0;i<*m;i++)
    {
      t = (double)i/(double)(*m);
      Au = inv_A_Pickands(*n, S, T, t);
      Ac = Au - (1.0 - t) * (invA0 - 1.0) - t * (invA1 - 1.0); // correction
      du = 1 / Au - Atheta[i];
      // dAinv =  Ac - 1 / Atheta[i];
      dc = 1 / Ac - Atheta[i];
      stat[0] += dc * dc;
      stat[1] += du * du;
      // stat[1] += dAinv * dAinv;
    }
  stat[0] = stat[0] * (double)(*n)/(double)(*m);
  stat[1] = stat[1] * (double)(*n)/(double)(*m);
}

/**
 * Cramer-von Mises test statistic based on the CFG estimator
 * used in the GOF for EV copulas
 * stat = n/m \sum_{i=0}^{m-1} [An(i/m) - A_{theta_n}(i/m)]^2
 *
 * @param n sample size
 * @param m grid size
 * @param S unit Fréchet pseudo-observations
 * @param T unit Fréchet pseudo-observations
 * @param Atheta values of the fitted A at the grid points
 * @param stat value of the test statistic
 * @author Ivan Kojadinovic
 */
void cramer_vonMises_CFG(int *n, int *m, double *S,
			 double *T, double *Atheta,
			 double *stat)
{
  int i;
  double t, Au, Ac, dc, du,
    logA0 = log_A_CFG(*n, S, T, 0.0),
    logA1 = log_A_CFG(*n, S, T, 1.0);

  stat[0] = 0.0; stat[1] = 0.0;
  for (i=0;i<*m;i++)
    {
      t = (double) i / (double) (*m);
      Au = log_A_CFG(*n, S, T, t);
      Ac = Au - (1.0 - t) * logA0 - t * logA1; // endpoint corrected
      dc = exp(Ac) - Atheta[i];
      du = exp(Au) - Atheta[i];
      // dlogA = Ac - log(Atheta[i]);
      stat[0] += dc * dc;
      stat[1] += du * du;
      // stat[1] += dlogA * dlogA;
    }
  stat[0] = stat[0] * (double)(*n)/(double)(*m);
  stat[1] = stat[1] * (double)(*n)/(double)(*m);
}


/* wrapper */
/**
 * Wrapper for the Cramer-von Mises test statistics
 * used in the GOF for EV copulas
 * stat = n/m \sum_{i=0}^{m-1} [An(i/m) - A_{theta_n}(i/m)]^2
 *
 * @param n sample size
 * @param m grid size
 * @param S unit Fréchet pseudo-observations
 * @param T unit Fréchet pseudo-observations
 * @param Atheta values of the fitted A at the grid points
 * @param stat value of the test statistic
 * @param CFG if > 0 then CFG, else Pickands
 * @author Ivan Kojadinovic
 */
void cramer_vonMises_Afun(int *n, int *m, double *S,
			  double *T, double *Atheta,
			  double *stat, int *CFG)
{
  if (*CFG)
    cramer_vonMises_CFG(n, m, S, T, Atheta, stat);
  else
    cramer_vonMises_Pickands(n, m, S, T, Atheta, stat);
}
