
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

/* fn_argdef_phase: define all accepted arguments for the 'phase' function.
 */
static const fn_args fn_argdef_phase[] = {
  { "ph0",     FN_ARGTYPE_FLOAT, "0.0" },
  { "ph1",     FN_ARGTYPE_FLOAT, "0.0" },
  { "pivot",   FN_ARGTYPE_FLOAT, "0.5" },
  { "ppm",     FN_ARGTYPE_BOOL,  "0" },
  { "hz",      FN_ARGTYPE_BOOL,  "0" },
  { "inverse", FN_ARGTYPE_BOOL,  "0" },
  { NULL, '\0', NULL }
};

/* fn_execute_phase(): applies a single-dimension phase correction operation.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application.
 * @args: function argument string.
 */
int fn_execute_phase (datum *D, const int dim, const char *argstr) {
  /* declare variables to hold argument values. */
  real ph0, ph1, piv;
  int inv, ppm, hz;

  /* declare required variables:
   * @ph: temporary array of phase correction values.
   * @szk: size of phase correction array dimension.
   */
  hx_array ph;
  int szk;

  /* parse the function argument string. */
  if (!fn_scan_args(argstr, fn_argdef_phase,
         &ph0, &ph1, &piv, &ppm, &hz, &inv))
    throw("failed to parse phase arguments");

  /* check the dimension index. */
  if (dim < 0 || dim >= D->nd)
    throw("dimension index %d out of bounds [0,%u)", dim, D->nd);

  /* check that the dimension is frequency-domain. */
  if (!D->dims[dim].ft)
    throw("dimension %d is not frequency-domain", dim);

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
  szk = D->array.sz[dim];
  if (!hx_array_alloc(&ph, D->array.d, 1, &szk))
    throw("failed to allocate phasor array");

  /* compute the phase correction values. */
  if (!hx_array_phasor(&ph, dim, ph0, ph1, piv))
    throw("failed to compute phasor array");

  /* perform the phase correction operation. */
  if (!hx_array_mul_vector(&D->array, &ph, dim, &D->array))
    throw("failed to execute phase correction");

  /* free the temporary arrays. */
  hx_array_free(&ph);

  /* return success. */
  return 1;
}

