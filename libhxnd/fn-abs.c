
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

/* fn_execute_abs(): compute the magnitude of the array of a datum structure.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application.
 * @args: function argument string.
 */
int fn_execute_abs (datum *D, const int dim, const char *argstr) {
  /* declare a required variable:
   * @d: local copy of the dimension index.
   */
  unsigned int d;

  /* check that no dimension was specified. */
  if (dim >= 0)
    throw("dimension index specification not supported");

  /* compute the absolute value of the datum array. */
  if (!hx_array_norm(&D->array))
    throw("failed to compute absolute value");

  /* drop the imaginary datum array coefficients. */
  if (!hx_array_real(&D->array, DATUM_DIM_INVALID))
    throw("failed to drop imaginaries");

  /* loop over the dimensions, converting their status to real. */
  for (d = 0; d < D->nd; d++) {
    /* convert the dimension to real. */
    D->dims[d].cx = 0;

    /* invalidate the array algebraic dimension. */
    D->dims[d].d = DATUM_DIM_INVALID;
  }

  /* return success. */
  return 1;
}

