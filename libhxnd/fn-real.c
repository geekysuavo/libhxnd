
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

/* fn_real(): drop imaginaries from the array of a datum structure.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application, or -1.
 * @args: function argument definition array.
 */
int fn_real (datum *D, const int dim, const fn_arg *args) {
  /* declare a required variable:
   * @drm: array algebraic dimension to be removed.
   * @d: dimension loop counter.
   */
  unsigned int d, drm;

  /* base function behaviour on the dimension index. */
  if (dim < 0) {
    /* drop all imaginary components from the array. */
    if (!hx_array_real(&D->array, DATUM_DIM_INVALID))
      throw("failed to drop imaginaries");

    /* invalidate the algebraic dimension indices. */
    for (d = 0; d < D->nd; d++) {
      D->dims[d].d = DATUM_DIM_INVALID;
      D->dims[d].cx = 0;
    }
  }
  else if (dim < D->nd) {
    /* check that the dimension is not already real. */
    if (!D->dims[dim].cx || D->dims[dim].d == DATUM_DIM_INVALID)
      return 1;

    /* store the array dimension to be removed. */
    drm = D->dims[dim].d;

    /* drop the specified imaginary component from the array. */
    if (!hx_array_real(&D->array, drm))
      throw("failed to drop imaginaries in dimension %d", dim);

    /* loop over the datum dimensions. */
    for (d = 0; d < D->nd; d++) {
      /* skip dimensions that are already real. */
      if (D->dims[d].d == DATUM_DIM_INVALID)
        continue;

      /* check if the dimension metadata needs adjustment. */
      if (D->dims[d].d == drm) {
        /* invalidate the specified dimension index. */
        D->dims[d].d = DATUM_DIM_INVALID;
        D->dims[d].cx = 0;
      }
      else if (D->dims[d].d > drm) {
        /* reduce the array algebraic dimension index. */
        D->dims[d].d--;
      }
    }
  }
  else
    throw("dimension index %d out of bounds [0,%u)", dim, D->nd);

  /* return success. */
  return 1;
}

