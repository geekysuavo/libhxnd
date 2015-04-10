
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

/* fn_symm_cb(): vector callback function used to symmetrize each pair of
 * vectors along two given dimensions of an array.
 */
int fn_symm_cb (hx_array *x, hx_array *y, int *arr, int idx, va_list *vl) {
  /* declare a few required variables:
   */
  int i, j, k, idxu, idxl, idxmin, idxmax, n;
  real unorm, lnorm, f;

  /* compute the number of scalar elements along each dimension. */
  n = y->sz[0];

  /* loop over the first dimension of the array. */
  for (i = 0; i < n; i++) {
    /* loop over the 'upper triangle' of the array. */
    for (j = i + 1; j < n; j++) {
      /* compute the upper and lower triangular linear indices. */
      idxu = i + n * j;
      idxl = j + n * i;

      /* perform different comparisons based on the array's complex state. */
      if (hx_array_is_real(y)) {
        /* use the real intensities of the two elements. */
        unorm = y->x[idxu];
        lnorm = y->x[idxl];
      }
      else {
        /* compute the absolute values of the two elements. */
        unorm = hx_data_real_norm(y->x + idxu * y->n, y->n);
        lnorm = hx_data_real_norm(y->x + idxl * y->n, y->n);
      }

      /* determine which point has the smaller absolute value. */
      idxmin = (unorm < lnorm ? idxu : idxl) * y->n;
      idxmax = (unorm < lnorm ? idxl : idxu) * y->n;

      /* compute the element scaling factor. */
      f = (unorm < lnorm ? unorm / lnorm : lnorm / unorm);

      /* perform different replacements based on whether or not
       * the input array is real.
       */
      if (hx_array_is_real(y)) {
        /* perform direct replacement. */
        y->x[idxmax] = y->x[idxmin];
      }
      else {
        /* perform signed replacement. */
        for (k = 0; k < y->n; k++)
          y->x[idxmax + k] *= f;
      }
    }
  }

  /* return success. */
  return 1;
}

/* fn_symm(): symmetrize the array of a datum structure along two dimensions.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application, or -1.
 * @args: function argument definition array.
 */
int fn_symm (datum *D, const int dim, const fn_arg *args) {
  /* declare variables to hold argument values:
   * @dims: array of dimension indices to apply symmetrization to.
   * @ndims: size of the dimension index array. should be two.
   */
  int ndims, *dims = NULL;

  /* declare a few required variables:
   * @k1, @k2: array topological dimension indices for symmetrization.
   */
  int k1, k2;

  /* check that no dimension was specified. */
  if (dim >= 0)
    throw("single dimension index specification not supported");

  /* get the argument values from the argdef. */
  if (!fn_args_get_all(args, &dims) ||
      !fn_args_get_sizes(args, &ndims))
    throw("failed to get symmetrize arguments");

  /* check if no dimensions were specified. */
  if (!dims) {
    /* allocate a default integer array. */
    ndims = 2;
    dims = hx_index_alloc(ndims);
    dims[0] = 1;
    dims[1] = 2;
  }

  /* check that the correct number of dimensions was specified. */
  if (ndims != 2)
    throw("unsupported symmetrization dimension count (%d != 2)", ndims);

  /* subtract down the dimension indices. */
  dims[0]--;
  dims[1]--;

  /* check that the first dimension index is in bounds. */
  if (dims[0] < 0 || dims[0] >= D->nd)
    throw("first dimension index %d out of bounds [0,%u)", dims[0], D->nd);

  /* check that the first dimension index is in bounds. */
  if (dims[1] < 0 || dims[1] >= D->nd)
    throw("second dimension index %d out of bounds [0,%u)", dims[1], D->nd);

  /* get the array topological indices. */
  k1 = D->dims[dims[0]].k;
  k2 = D->dims[dims[1]].k;

  /* check that the dimensions have the same point count. */
  if (D->array.sz[k1] != D->array.sz[k2])
    throw("symmetrization requires square planes");

  /* symmetrize each submatrix of the array. */
  if (!hx_array_foreach_matrix(&D->array, k1, k2, &fn_symm_cb))
    throw("failed to apply symmetrization");

  /* free the index array. */
  hx_index_free(dims);

  /* return success. */
  return 1;
}

