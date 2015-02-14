
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

/* fn_mirror_cb(): vector callback function used to mirror each vector along
 * a given dimension of an array.
 */
int fn_mirror_cb (hx_array *x, hx_array *y, int *arr, int idx, va_list *vl) {
  /* declare a few required variables:
   * @i: array coefficient index of the lower off-zero point.
   * @j: array coefficient index of the upper off-zero point.
   * @n: number of scalar elements in the sliced array.
   * @inorm, @jnorm: absolute values of the scalars at @i and @j.
   * @idxmin: array coefficient index of the smaller scalar.
   * @idxmax: array coefficient index of the larger scalar.
   * @f: scaling factor based on the absolute values.
   */
  int i, j, k, idxmin, idxmax, n;
  real inorm, jnorm, f;

  /* compute the number of scalar elements in the array. */
  n = y->len / y->n;

  /* initialize the starting indices. */
  if (n % 2) {
    /* odd number of elements. */
    i = n / 2 - 1;
    j = n / 2 + 1;
  }
  else {
    /* even number of elements. */
    i = n / 2 - 1;
    j = n / 2;
  }

  /* loop over the off-zero elements of the array. */
  while (i >= 0 && j < n) {
    /* perform different comparisons based on the array's complex state. */
    if (hx_array_is_real(y)) {
      /* use the real intensities of the two elements. */
      inorm = y->x[i];
      jnorm = y->x[j];
    }
    else {
      /* compute the absolute values of the two elements. */
      inorm = hx_data_real_norm(y->x + i * y->n, y->d, y->n);
      jnorm = hx_data_real_norm(y->x + j * y->n, y->d, y->n);
    }

    /* determine which point has the smaller absolute value. */
    idxmin = (inorm < jnorm ? i : j) * y->n;
    idxmax = (inorm < jnorm ? j : i) * y->n;

    /* compute the element scaling factor. */
    f = (inorm < jnorm ? inorm / jnorm : jnorm / inorm);

    /* perform different replacements based on whether
     * or not the input array is real.
     */
    if (hx_array_is_real(y)) {
      /* perform direct replacement. */
      y->x[idxmax] = y->x[idxmin];
    }
    else {
      /* perform signed replacement. effectively, scaling the value by
       * the ratio of the norms replaces its value with the point
       * with the smaller absolute value, but retains the sign
       * of the original value.
       */
      for (k = 0; k < y->n; k++)
        y->x[idxmax + k] *= f;
    }

    /* move the indices to the next pair of elements. */
    i--;
    j++;
  }

  /* return success. */
  return 1;
}

/* fn_mirror(): mirror the array of a datum structure along a given dimension.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application, or -1.
 * @args: function argument definition array.
 */
int fn_mirror (datum *D, const int dim, const fn_arg *args) {
  /* check the dimension index. */
  if (dim < 0 || dim >= D->nd)
    throw("dimension index %d out of bounds [0,%u)", dim, D->nd);

  /* execute the mirroring vector operation. */
  if (!hx_array_vector_op(&D->array, D->dims[dim].k, &fn_mirror_cb))
    throw("failed to perform mirroring operation");

  /* return success. */
  return 1;
}

