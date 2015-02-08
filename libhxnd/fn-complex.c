
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

/* fn_complex(): add imaginaries into the array of a datum structure.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application, or -1.
 * @args: function argument definition array.
 */
int fn_complex (datum *D, const int dim, const fn_arg *args) {
  /* declare a required variable:
   * @allcx: flag that it set if all dimensions are complex.
   * @d: dimension loop counter.
   */
  unsigned int d, dnew;
  int allcx;

  /* loop over the datum dimensions. */
  for (d = 0, allcx = 1; d < D->nd; d++) {
    /* check if the current dimension is real. */
    if (!D->dims[d].cx || D->dims[d].d == DATUM_DIM_INVALID)
      allcx = 0;
  }

  /* if all dimensions are already complex, return success. */
  if (allcx)
    return 1;

  /* base function behaviour on the dimension index. */
  if (dim < 0) {
    /* store the maximum algebraic dimension index currently available. */
    dnew = D->array.d - 1;

    /* make every array dimension complex. */
    if (!hx_array_resize(&D->array, D->array.k, D->array.k, D->array.sz))
      throw("failed to complex promote all array dimensions");

    /* loop over the datum dimensions. */
    for (d = 0; d < D->nd; d++) {
      /* check if the current dimension was real. */
      if (!D->dims[d].cx) {
        /* complex-promote the dimension. */
        D->dims[d].cx = 1;
        D->dims[d].d = ++dnew;
      }
    }
  }
  else if (dim < D->nd) {
    /* return success if the dimension is already complex. */
    if (D->dims[dim].cx && D->dims[dim].d > DATUM_DIM_INVALID)
      return 1;

    /* store the maximum algebraic dimension index currently available. */
    dnew = D->array.d - 1;

    /* make the specified array dimension complex. */
    if (!hx_array_resize(&D->array, D->array.d + 1, D->array.k, D->array.sz))
      throw("failed to complex promote array dimension");

    /* complex-promote the dimension. */
    D->dims[dim].cx = 1;
    D->dims[dim].d = ++dnew;
  }
  else
    throw("dimension index %d out of bounds [0,%u)", dim, D->nd);

  /* return success. */
  return 1;
}

