
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

/* include the processing function header. */
#include <hxnd/fn.h>
#include <hxnd/fn-handlers.h>

/* fn_baseline_cb(): vector callback function used to baseline-correct
 * each vector along a given dimension of an array.
 */
int fn_baseline_cb (hx_array *x, hx_array *y,
                    hx_index idx, int pidx,
                    va_list *vl) {
  /* extract the varargs. */
  hx_array *w = va_arg(*vl, hx_array*);
  hx_array *y0 = va_arg(*vl, hx_array*);
  real smooth = (real) va_arg(*vl, double);

  /* compute the baseline weights. */
  if (!hx_baseline_weight(y, w))
    throw("failed to compute baseline weights");

  /* solve for the smoothed baseline. */
  if (!hx_baseline(y, w, smooth, y0))
    throw("failed to compute smoothed baseline");

  /* compute the baseline-corrected trace vector. */
  if (!hx_array_add_array(y, y0, -1.0, y))
    throw("failed to subtract baseline");

  /* return success. */
  return 1;
}

/* fn_baseline(): applies a first-dimension baseline correction operation.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application, or -1.
 * @args: function argument definition array.
 */
int fn_baseline (datum *D, const int dim, const fn_arg *args) {
  /* declare variables to hold argument values. */
  real smooth; /* ;) */

  /* declare required variables:
   * @w: array holding each weight vector.
   * @b: array holding each baseline vector.
   * @d: array algebraic dimensionality.
   * @k: array topological dimension index.
   * @szk: array baseline/weight vector length.
   */
  hx_array w, b;
  int d, k, szk;

  /* get the argument values from the argdef array */
  if (!fn_args_get_all(args, &smooth))
    throw("failed to get baseline arguments");

  /* check that no dimension was specified. */
  if (dim >= 0)
    throw("dimension index specification not supported");

  /* check that the dimension is frequency-domain. */
  if (!D->dims[0].ft)
    throw("first dimension is not frequency-domain");

  /* get the dimensionality of the array. */
  d = D->array.d;

  /* get the dimension index and size of the correction. */
  k = D->dims[0].k;
  szk = D->array.sz[k];

  /* allocate a vector of weights. */
  if (!hx_array_alloc(&w, 0, 1, &szk))
    throw("failed to allocate weight (0, 1)-array");

  /* allocate a vector of baseline values. */
  if (!hx_array_alloc(&b, d, 1, &szk))
    throw("failed to allocate baseline (%d, 1)-array", d);

  /* apply the required baseline corrections. */
  if (!hx_array_foreach_vector(&D->array, k, &fn_baseline_cb,
                               &w, &b, smooth))
    throw("failed to apply baseline correction");

  /* free the temporary vectors. */
  hx_array_free(&w);
  hx_array_free(&b);

  /* return success. */
  return 1;
}

