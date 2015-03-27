
/* hxnd: A framework for n-dimensional hypercomplex calculations for NMR.
 * Copyright (C) 2014  Bradley Worley  <geekysuavo@gmail.com>.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to:
 *
 *   Free Software Foundation, Inc.
 *   51 Franklin Street, Fifth Floor
 *   Boston, MA  02110-1301, USA.
 */

/* include the baseline header. */
#include <hxnd/hx-baseline.h>

/* hx_baseline_weight(): build an array of weights based on a rough guess
 * from a vector of hypercomplex spectral intensities.
 * @x: pointer to the array of spectral intensities.
 * @w: pointer to the array to store weights into.
 */
int hx_baseline_weight (hx_array *x, hx_array *w) {
  /* declare required variables:
   * @mu: mean value of the y-vector.
   * @sigma: standard deviation of the y-vector.
   * @th: current weight selection threshold.
   * @y: array of first-derivative magnitudes.
   */
  real mu, sigma, th;
  int n, N, W, Wprev;
  hx_array y;

  /* compute the length of the input vector. */
  N = x->len / x->n;

  /* allocate a temporary norm array. */
  if (!hx_array_copy(&y, x))
    throw("failed to allocate norm array");

  /* compute the first differences of the input array. */
  for (n = N - 1; n >= 1; n--)
    hx_data_add(y.x + y.n * n,
                y.x + y.n * (n - 1),
                y.x + y.n * n,
                -1.0, x->d, x->n);

  /* zero the first array element. */
  hx_data_zero(y.x, y.n);

  /* compute the norm of each vector element. */
  if (!hx_array_norm(&y) ||
      !hx_array_resize(&y, 0, y.k, y.sz))
    throw("failed to compute norm array");

  /* square each element of the norm vector. */
  for (n = 0; n < N; n++)
    y.x[n] *= y.x[n];

  /* fill the weight vector with ones. */
  W = N;
  if (!hx_array_fill(w, 1.0))
    throw("failed to initialize weights");

  /* loop until we have estimated the baseline points. */
  do {
    /* store the previous weight count. */
    Wprev = W;

    /* compute the mean of the weighed data points. */
    for (n = 0, mu = 0.0; n < N; n++)
      mu += w->x[n] * y.x[n];

    /* normalize the mean. */
    mu /= (real) Wprev;

    /* compute the standard deviation of the weighted data points. */
    for (n = 0, sigma = 0.0; n < N; n++)
      sigma += w->x[n] * pow(y.x[n] - mu, 2.0);

    /* normalize the standard deviation. */
    sigma = sqrt(sigma / (real) (Wprev - 1));

    /* compute the threshold intensity. */
    th = mu + 2.0 * sigma;

    /* recompute the weights. */
    for (n = W = 0; n < N; n++) {
      /* check if the current intensity is below the threshold. */
      if (y.x[n] <= th) {
        /* store the weight. */
        w->x[n] = 1.0;
        W++;
      }
      else {
        /* store the weight. */
        w->x[n] = 0.0;
      }
    }
  } while (W != Wprev);

  /* free the temporary vector. */
  hx_array_free(&y);

  /* loop over the weights to erode the values. */
  for (n = 1; n < N - 1; n++) {
    /* remove lone positive weights. */
    if (w->x[n - 1] == 0.0 && w->x[n + 1] == 0.0)
      w->x[n] = 0.0;

    /* remove lone negative weights. */
    if (w->x[n - 1] == 1.0 && w->x[n + 1] == 1.0)
      w->x[n] = 1.0;
  }

  /* return success. */
  return 1;
}

/* hx_baseline(): compute the whittaker-smoothed baseline vector from a
 * vector of hypercomplex spectral intensities and its corresponding
 * vector of weights, along with a smoothness factor.
 *
 * this function implements a niche sparse cholesky decomposition to
 * compute the final vector of baseline intensities. the known structure
 * of the sparse matrix (symmetric, positive definite, tridiagonal) is
 * taken advantage of.
 *
 * @x: pointer to the array of spectral intensities.
 * @w: pointer to an array weights (0: interpolate <-> 1: fit)
 * @lambda: smoothing factor. larger values create smoother fits.
 * @x0: pointer to the array to store the fitted baseline into.
 */
int hx_baseline (hx_array *x, hx_array *w, real lambda, hx_array *x0) {
  /* declare required variables:
   * @a: main diagonal of the regression design matrix.
   * @b: first subdiagonal of the design matrix.
   * @n: vector element loop counter.
   * @N: vector element count.
   */
  real lambda0, *xcur, *xpre, *xnxt;
  hx_array a, b;
  int n, N;

  /* compute the length of the vector to smooth. */
  N = x->len / x->n;

  /* allocate the regression arrays. */
  if (!hx_array_alloc(&a, 0, 1, &N) ||
      !hx_array_alloc(&b, 0, 1, &N))
    throw("failed to allocate design matrix");

  /* compute the base smoothing factor. */
  for (n = 0, lambda0 = 0.0; n < N; n++)
    lambda0 += (real) w->x[n];

  /* compute the elements of the regression arrays. */
  for (n = 0; n < N; n++) {
    /* store the diagonal element. */
    if (n == 0 || n == N - 1)
      a.x[n] = w->x[n] + 1.0 * lambda0 * lambda;
    else
      a.x[n] = w->x[n] + 2.0 * lambda0 * lambda;

    /* store the subdiagonal element. */
    if (n > 0)
      b.x[n] = -1.0 * lambda0 * lambda;
  }

  /* replace the regression array with its cholesky factors. */
  for (n = 0; n < N; n++) {
    /* compute the diagonal element. */
    a.x[n] = sqrt(a.x[n] - pow(b.x[n], 2.0));

    /* compute the subdiagonal element. */
    if (n < N - 1)
      b.x[n + 1] /= a.x[n];
  }

  /* perform the first half of the sparse cholesky solution. */
  for (n = 0; n < N; n++) {
    /* compute the current and previous data point addresses. */
    xcur = x0->x + x->n * n;
    xpre = x0->x + x->n * (n - 1);

    /* copy the data point from the input vector. */
    hx_data_add(NULL, x->x + n * x->n, xcur,
                w->x[n], x->d, x->n);

    /* scale the data point by the diagonal. */
    hx_data_add(NULL, xcur, xcur,
                1.0 / a.x[n],
                x->d, x->n);

    /* subtract the off-diagonal when possible. */
    if (n > 0)
      hx_data_add(xcur, xpre, xcur,
                  -b.x[n] / a.x[n],
                  x->d, x->n);
  }

  /* perform the second half of the sparse cholesky solution. */
  for (n = N - 1; n >= 0; n--) {
    /* compute the current and next data point addresses. */
    xcur = x0->x + x->n * n;
    xnxt = x0->x + x->n * (n + 1);

    /* scale the data point by the diagonal. */
    hx_data_add(NULL, xcur, xcur,
                1.0 / a.x[n],
                x->d, x->n);

    /* subtract the off-diagonal when possible. */
    if (n < N - 1)
      hx_data_add(xcur, xnxt, xcur,
                  -b.x[n + 1] / a.x[n],
                  x->d, x->n);
  }

  /* free the regression arrays. */
  hx_array_free(&a);
  hx_array_free(&b);

  /* return success. */
  return 1;
}

