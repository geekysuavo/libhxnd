
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

/* fn_tilt_cb(): vector callback function used to shift each vector
 * along a given dimension of an array, relative to another dimension.
 */
int fn_tilt_cb (hx_array *x, hx_array *y, int *arr, int idx, va_list *vl) {
  /* declare a few required variables:
   * @f: tilt factor to apply.
   * @n: number of points to shift.
   * @k: dimension to use as a reference.
   * @d: algebraic basis index for shifting.
   */
  real f, n;
  int d, k;

  /* extract the varargs. */
  f = (real) va_arg(*vl, double);
  d = va_arg(*vl, int);
  k = va_arg(*vl, int);

  /* compute the shift amount. */
  n = f * ((real) x->sz[k] / 2.0 - (real) arr[k]);

  /* shift the sliced vector. */
  if (!hx_array_fshift(y, d, 0, n))
    throw("failed to execute fractional shift");

  /* return success. */
  return 1;
}

/* fn_tilt(): shift the vectors of the array of a datum structure in order to
 * achieve a given 'tilt angle'.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application, or -1.
 * @args: function argument definition array.
 */
int fn_tilt (datum *D, const int dim, const fn_arg *args) {
  /* declare variables to hold argument values:
   * @dims: array of dimension indices to apply tilting to.
   * @ndims: size of the dimension index array. should be two.
   * @angle: angle of tilt to apply to the dimension pair.
   */
  int ndims, *dims = NULL;
  real angle;

  /* declare a few required variables:
   * @d1: array algebraic dimension index to shift.
   * @k1: array topological dimension index to shift.
   * @k2: array topological dimension index for computing shift amount.
   */
  int d1, k1, k2;

  /* check that no dimension was specified. */
  if (dim >= 0)
    throw("single dimension index specification not supported");

  /* get the argument values from the argdef. */
  if (!fn_args_get_all(args, &angle, &dims) ||
      !fn_args_get_sizes(args, NULL, &ndims))
    throw("failed to get tilt arguments");

  /* check if no dimensions were specified. */
  if (!dims) {
    /* allocate a default integer array. */
    ndims = 2;
    dims = hx_array_index_alloc(ndims);
    dims[0] = 1;
    dims[1] = 2;
  }

  /* check that the correct number of dimensions was specified. */
  if (ndims != 2)
    throw("unsupported tilt dimension count (%d != 2)", ndims);

  /* subtract down the dimension indices. */
  dims[0]--;
  dims[1]--;

  /* check that the first dimension index is in bounds. */
  if (dims[0] < 0 || dims[0] >= D->nd)
    throw("first dimension index %d out of bounds [0,%u)", dims[0], D->nd);

  /* check that the first dimension index is in bounds. */
  if (dims[1] < 0 || dims[1] >= D->nd)
    throw("second dimension index %d out of bounds [0,%u)", dims[1], D->nd);

  /* check that the shift dimension is complex. */
  if (D->dims[dims[0]].d == DATUM_DIM_INVALID)
    throw("shift dimension must be complex for tilt");

  /* if unspecified, automatically compute the tilt angle. */
  if (angle == 0.0)
    angle = D->dims[dims[1]].width / D->dims[dims[0]].width;

  /* set the shift and count dimension index. */
  d1 = D->dims[dims[0]].d;
  k1 = D->dims[dims[0]].k;
  k2 = D->dims[dims[1]].k;

  /* check that the tilt angle is either specified or valid. */
  if (angle == 0.0)
    throw("invalid or unspecified tilt angle");

  /* scale the tilt angle into point units. */
  angle *= (real) D->dims[dims[0]].sz / (real) D->dims[dims[1]].sz;

  /* apply the required shift operations. */
  if (!hx_array_foreach_vector(&D->array, k1, &fn_tilt_cb, angle, d1, k2))
    throw("failed to apply tilt operation");

  /* free the index array. */
  free(dims);

  /* return success. */
  return 1;
}

