
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

/* fn_shift(): applies a single-dimension point-shift operation.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application, or -1.
 * @args: function argument definition array.
 */
int fn_shift (datum *D, const int dim, const fn_arg *args) {
  /* declare variables to hold argument values.
   */
  int pts, sec, ppm, hz, doround, nset;

  /* declare required variables:
   * @famt: floating-point shift amount.
   * @iamt: integer shift amount.
   */
  real famt, fsz;
  int iamt;

  /* get the argument values from the argdef array. */
  if (!fn_args_get_all(args, &pts, &sec, &ppm, &hz, &doround, &famt))
    throw("failed to get shift arguments");

  /* check the dimension index. */
  if (dim < 0 || dim >= D->nd)
    throw("dimension index %d out of bounds [0,%u)", dim, D->nd);

  /* determine how many unit options were set. */
  nset = pts + sec + ppm + hz;

  /* if no options were set, default to points. */
  if (nset == 0)
    pts = 1;

  /* if multiple options were set, throw an exception. */
  if (nset > 1)
    throw("multiple unit options set");

  /* check that frequency units were used on frequency-domain data. */
  if (D->dims[dim].ft && sec)
    throw("cannot specify seconds on frequency-domain dimensions");

  /* check that time units were used on time-domain data. */
  if (!D->dims[dim].ft && (hz || ppm))
    throw("cannot specify hz/ppm on time-domain dimensions");

  /* store the dimension size in floating point. */
  fsz = (real) D->dims[dim].sz;

  /* compute the shift amount. */
  if (sec) {
    /* convert the value in seconds to points. */
    famt *= D->dims[dim].width;
  }
  else if (hz) {
    /* convert the value in hertz to points. */
    famt *= (fsz / D->dims[dim].width);
  }
  else if (ppm) {
    /* convert the value in ppm to points. */
    famt *= (D->dims[dim].carrier * fsz / D->dims[dim].width);
  }

  /* check if fractional shifts are allowed. */
  if (doround) {
    /* no. compute the integer-valued shift amount. */
    iamt = (int) famt;

    /* perform the integer shift. */
    if (!hx_array_shift(&D->array, D->dims[dim].k, iamt))
      throw("failed to perform integer shift");
  }
  else {
    /* check that the shift dimension is complex. */
    if (D->dims[dim].d == DATUM_DIM_INVALID)
      throw("fractional shift dimension must be complex");

    /* perform the fractional shift. */
    if (!hx_array_fshift(&D->array, D->dims[dim].d, D->dims[dim].k, famt))
      throw("failed to perform fractional shift");
  }

  /* return success. */
  return 1;
}

