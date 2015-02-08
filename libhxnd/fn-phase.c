
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

/* fn_phase(): applies a single-dimension phase correction operation.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application, or -1.
 * @args: function argument definition array.
 */
int fn_phase (datum *D, const int dim, const fn_arg *args) {
  /* declare variables to hold argument values. */
  real ph0, ph1, piv;
  int inv, ppm, hz;

  /* declare required variables:
   * @ph: temporary array of phase correction values.
   * @szk: size of phase correction array dimension.
   */
  hx_array ph;
  int k, szk;

  /* get the argument values from the argdef array */
  if (!fn_args_get_all(args, &ph0, &ph1, &piv, &ppm, &hz, &inv))
    throw("failed to get phase arguments");

  /* check the dimension index. */
  if (dim < 0 || dim >= D->nd)
    throw("dimension index %d out of bounds [0,%u)", dim, D->nd);

  /* check that the dimension is frequency-domain. */
  if (!D->dims[dim].ft)
    throw("dimension %d is not frequency-domain", dim);

  /* check that the dimension is complex. */
  if (!D->dims[dim].cx)
    throw("dimension %d is not complex", dim);

  /* check that no more than one unit was specified. */
  if (ppm && hz)
    throw("multiple unit options set");

  /* convert the pivot from ppm to hertz, if required. */
  if (ppm)
    piv = piv * D->dims[dim].carrier;

  /* convert the pivot from hertz to fractional index, if required. */
  if (ppm || hz)
    piv = ((piv - D->dims[dim].offset) / D->dims[dim].width) + 0.5;

  /* check that the pivot is within bounds. */
  if (piv < 0.0 || piv > 1.0)
    throw("pivot value %.3f is out of bounds [0,1]", piv);

  /* invert the phase correction values, if requested. */
  if (inv) {
    /* yep. negate @ph0 and @ph1. */
    ph0 *= -1.0;
    ph1 *= -1.0;
  }

  /* convert the phase correction values to radians. */
  ph0 *= (M_PI / 180.0);
  ph1 *= (M_PI / 180.0);

  /* allocate a temporary phase correction vector. */
  k = D->dims[dim].k;
  szk = D->array.sz[k];
  if (!hx_array_alloc(&ph, D->array.d, 1, &szk))
    throw("failed to allocate phasor array");

  /* compute the phase correction values. */
  if (!hx_array_phasor(&ph, D->dims[dim].d, ph0, ph1, piv))
    throw("failed to compute phasor array");

  /* perform the phase correction operation. */
  if (!hx_array_mul_vector(&D->array, &ph, k, &D->array))
    throw("failed to execute phase correction");

  /* free the temporary arrays. */
  hx_array_free(&ph);

  /* return success. */
  return 1;
}

