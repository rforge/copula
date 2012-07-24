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
 * @file   ecIndepTest.c
 * @author Ivan Kojadinovic
 * @date   December 2007
 *
 * @brief  Multivariate independence test based on the empirical
 *         copula process as proposed by Christian Genest and Bruno
 *         R�millard (2004), Test 13:2, pages 335-369.
 *
 */

#include <R.h>
#include <Rmath.h>
#include "set_utils.h"
#include "indepTests.h"

/// Temporary array J
void J_u(int n, int p, const double R[], double *J)
{
  int i, j, l, m;

  m=0;
  for (j=0;j<p;j++)
    for (l=0;l<n;l++)
      for (i=0;i<n;i++)
	J[m++] = 1.0 - fmax2(R[i + n*j], R[l + n*j])/n;

}

/**
 * Simulate the distribution of TA (up to subsets of cardinality p)
 * and of the global statistic
 *
 * @param n sample size
 * @param N number of simulations
 * @param p dimension
 * @param m max cardinality of A
 * @param TA0 values of TA under independence (N repetitions)
 * @param G0 values of the global stat. under independence (N repetitions)
 * @param subset subsets of {1,...,p} in binary notation (int) whose card. is
 *               between 2 and m in "natural" order
 * @param subset_char similar, for printing
 * @param fisher0 p-values � la Fisher
 * @param tippett0 p-values � la Tippett
 * @param verbose display progress bar if > 0
 * @author Ivan Kojadinovic
 */
void simulate_empirical_copula(int *n, int *N, int *p, int *m, double *TA0,
			       double *G0, int *subset, char **subset_char,
			       double *fisher0, double *tippett0, int *verbose)
{
  int i, j, k, index, sb[1], count;
  double *R = Calloc((*n) * (*p), double);
  double *J = Calloc((*n) * (*n) * (*p), double);
  double *K = Calloc((*n) * (*p), double);
  double *L = Calloc(*p, double);
  double pvalue, r;

  /* number of subsets */
  *sb = (int)sum_binom(*p,*m);

  /* partial power set in "natural" order */
  k_power_set(p, m, subset);

  /* convert partial power set to char for printing */
  k_power_set_char(p, sb, subset, subset_char);

  /* Simulate the distribution of TA under independence */

  GetRNGstate();

  /* N repetitions */
  for (k=0;k<*N;k++)
    {
      /* generate data */
      for (j=0;j<*p;j++)
	{
	  for (i=0;i<*n;i++)
	    R[i + (*n) * j] = i+1;

	  /* permutation = random ranks in column j */
	  for (i=*n-1;i>=0;i--)
	    {
	      r = R[j * (*n) + i];
	      index = (int)((i + 1) *  unif_rand());
	      R[j * (*n) + i] = R[j * (*n) + index];
	      R[j * (*n) + index] = r;
	    }
	}

      /* compute arrays J, K, L */
      J_u(*n, *p, R, J);
      K_array(*n, *p, J, K);
      L_array(*n, *p, K, L);

      /* for subsets i of cardinality greater than 1 and lower than m */
      for (i=*p+1;i<*sb;i++)
	TA0[k + (*N) * (i - *p - 1)] =  M_A_n(*n, *p, J, K, L, subset[i]);

      /* global stat under independence*/
      G0[k] = I_n(*n, *p, J, K, L);

      if (*verbose)
	progressBar(k, *N, 70);
    }

  PutRNGstate();

  /* compute W � la Fisher and � la Tippett from TA0 */
  for (k=0;k<*N;k++)
    {
      fisher0[k] = 0.0;
      tippett0[k] = 1.0;
      for (i=0;i<*sb-*p-1;i++)
	{
	  /* p-value */
	  count = 0;
	  for (j=0;j<*N;j++)
	    if (TA0[j + (*N) * i] >= TA0[k + (*N) * i])
	      count ++;
	  pvalue = (double)(count + 0.5)/(*N + 1.0);
	  fisher0[k] -= 2*log(pvalue);
	  tippett0[k] = fmin2(tippett0[k],pvalue);
	}
    }

  Free(R);

  Free(J);
  Free(K);
  Free(L);
}

/*****************************************************************************

  Compute TA for subsets of cardinality 2 to p
  R: ranks
  n: sample size
  p: number of variables
  m: max. cardinality of subsets of {1,...,p}
  TA0: simulated values of TA under independence (size: N * (sb - p - 1))
  N: number of repetitions (nrows TA0, fisher0, tippett0)
  subset: subsets of {1,...,p} in binary notation (int) whose card. is
  between 2 and m in "natural" order
  TA: test statistics (size: sum.bin. - p - 1)
  pval: corresponding p-values (size = sum.bin. - p - 1)
  fisher: pvalue � la Fisher
  tippett: pvalue � la Tippett

******************************************************************************/

/**
 * Computes the statistcs TA (up to subsets of cardinality p)
 * Computes the global statistic In
 *
 * @param R multivariate ranks
 * @param n sample size
 * @param p dimension
 * @param m maximum cardinality of subsets of {1,...,p}
 * @param TA0 simulated values of TA under independence (size: N * (sb - p - 1))
 * @param G0 simulated values of In under independence
 * @param N number of repetitions (nrows TA0, fisher0, tippett0)
 * @param subset subsets of {1,...,p} in binary notation (int) whose card. is
 *               between 2 and m in "natural" order
 * @param TA test statistics (size: sum.bin. - p - 1)
 * @param G global statistic
 * @param pval p-values (size = sum.bin. - p - 1) for TA
 * @param fisher pvalue � la Fisher
 * @param tippett pvalue � la Tippett
 * @param globpval pvalue of In
 * @param fisher0 p-values � la Fisher under the null
 * @param tippett0 p-values � la Tippett under the null
 * @author Ivan Kojadinovic
 */
void empirical_copula_test(double *R, int *n, int *p, int *m, double *TA0, double *G0,
			   int *N, int *subset, double *TA, double *G, double *pval,
			   double *fisher, double *tippett, double *globpval,
			   double *fisher0, double *tippett0)
{
  int i, k, count, sb = (int)sum_binom(*p,*m);

  double *J = Calloc((*n) * (*n) * (*p), double);
  double *K = Calloc((*n) * (*p), double);
  double *L = Calloc(*p, double);

  /* compute arrays J, K, L */
  J_u(*n, *p, R, J);
  K_array(*n, *p, J, K);
  L_array(*n, *p, K, L);

  /* compute W from the current data */
  *fisher = 0.0;
  *tippett = 1.0;

  /* for subsets i of cardinality greater than 1 */
  for (i=0;i<sb-*p-1;i++)
    {
      TA[i] = M_A_n(*n, *p, J, K, L, subset[i]);

      /* p-value */
      count = 0;
      for (k=0;k<*N;k++)
	if (TA0[k + (*N) * i] >= TA[i])
	  count ++;
      pval[i] = (double)(count + 0.5)/(*N + 1.0);

      *fisher -= 2*log(pval[i]);
      *tippett = fmin2(*tippett,pval[i]);
    }

  /* p-values of the Fisher and Tippett statistics */
  count = 0;
  for (k=0;k<*N;k++)
    if (fisher0[k] >= *fisher)
      count ++;
  *fisher = (double)(count + 0.5)/(*N + 1.0);

  count = 0;
  for (k=0;k<*N;k++)
    if (tippett0[k] <= *tippett)
      count ++;
  *tippett = (double)(count + 0.5)/(*N + 1.0);


  /* compute global stat from the current data and the corresponding p-value*/
  *G = I_n(*n, *p, J, K, L);

  count = 0;
  for (k=0;k<*N;k++)
    if (G0[k] >= *G)
      count ++;
  *globpval = (double)(count + 0.5)/(*N + 1.0);

  Free(J);
  Free(K);
  Free(L);
}


