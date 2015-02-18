
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

/* include the filter and window function headers. */
#include <hxnd/hx-filter.h>
#include <hxnd/hx-window.h>

/* hx_sinc(): compute the value of the sin(pi*x)/(pi*x) function.
 * @x: argument of the sinc() function.
 */
real hx_sinc (real x) {
  /* return the value of the sinc() function. */
  return (x == 0.0 ? 1.0 : sin(M_PI * x) / (M_PI * x));
}

/* hx_filter_fir_alloc(): construct a real finite impulse response filter
 * coefficient vector based on filter length and transition frequencies.
 * @b: pointer to the array to store the coefficients into.
 * @M: filter order. one less than the coefficient count..
 * @ft: normalized filter transition frequency.
 * @inv: whether to construct a band-pass (0) or band-stop (1) filter.
 */
int hx_filter_fir_alloc (hx_array *b, int M, real ft, int inv) {
  /* declare a few required variables:
   * @bi: currently indexed filter coefficient, unscaled.
   * @x: fractional coefficient array index.
   * @i: integral coefficient array index.
   */
  real bi, x;
  int i;

  /* check that the filter order supports inversion. */
  if (inv && M % 2)
    throw("band-stop filter must have even order");

  /* check that the transition frequency is within bounds. */
  if (ft < 0 || ft > 0.5)
    throw("transition frequency %.4f out of bounds [0,0.5]", ft);

  /* initialize the filter coefficients with a blackman window. */
  if (!hx_window_black(b, 0, M + 1))
    throw("failed to initialize filter coefficients");

  /* loop over the filter coefficients. */
  for (i = 0; i <= M; i++) {
    /* compute the fractional index. */
    x = (real) (i - M / 2);

    /* compute the filter coefficient. */
    bi = (inv ? -1.0 : 1.0) * 2.0 * ft * hx_sinc(2.0 * ft * x);

    /* determine whether to construct a pass or a stop filter. */
    if (inv && x == 0.0)
      bi = 1.0 - 2.0 * ft;

    /* scale the coefficient into the array. */
    b->x[i] *= bi;
  }

  /* return success. */
  return 1;
}

/* hx_filter_firfn(): per-vector callback function for hx_filter_fir()
 */
int hx_filter_firfn (hx_array *x, hx_array *y,
                     int *arr, int idx,
                     va_list *vl) {
  /* declare a few required variables:
   * @i: output array coefficient index.
   * @m: filter coefficient index.
   * @M: filter order.
   */
  int i, j, m, M;

  /* extract the vararg. */
  hx_array *b = va_arg(*vl, hx_array*);

  /* get the filter order. */
  M = b->len - 1;

  /* loop backwards over the array values. */
  for (i = y->len - y->n - 1; i >= 0; i -= y->n) {
    /* scale the current term. */
    for (j = 0; j < y->n; j++)
      y->x[i + j] *= b->x[0];

    /* loop over the filter convolution sum. */
    for (m = 1; m <= M; m++) {
      /* skip terms hanging off the side of the array. */
      if (i - m * y->n < 0)
        continue;

      /* sum the current scaled term. */
      for (j = 0; j < y->n; j++)
        y->x[i + j] += y->x[i - m * y->n + j] * b->x[m];
    }
  }

  /* return success. */
  return 1;
}

/* hx_filter_fir(): apply a constructed set of finite impulse response filter
 * coefficients to a dimension of an array.
 * @x: pointer to the array to apply the filter to.
 * @k: topological dimension index to filter.
 * @b: array of filter coefficients.
 */
int hx_filter_fir (hx_array *x, int k, hx_array *b) {
  /* check that the dimension index is in bounds. */
  if (k < 0 || k >= x->k)
    throw("dimension index %d out of bounds [0,%d)", k, x->k);

  /* check that the coefficient array is correctly shaped. */
  if (b->d != 0 || b->k != 1 || b->sz[0] >= x->sz[k])
    throw("invalid fir filter coefficient array");

  /* apply the filter coefficients to each k-mode vector of the array. */
  if (!hx_array_foreach_vector(x, k, &hx_filter_firfn, b))
    throw("failed to apply fir filter coefficients");

  /* return success. */
  return 1;
}

