
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

/* fn_argdef_scale: define all accepted arguments for the 'scale' function.
 */
static const fn_args fn_argdef_scale[] = {
  { "first",   FN_ARGTYPE_FLOAT,  "1.0" },
  { "factor",  FN_ARGTYPE_FLOAT,  "1.0" },
  { "invert",  FN_ARGTYPE_BOOL,   "0" },
  { NULL, '\0', NULL }
};

/* fn_execute_scale(): scales a datum structure by a constant factor.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application.
 * @args: function argument string.
 */
int fn_execute_scale (datum *D, const int dim, const char *argstr) {
  /* declare variables to hold argument values.
   * @hxscale: hypercomplex scalar for scaling..
   * @f0: first point scaling factor.
   * @fscale: scaling factor.
   * @atmp: temporary array.
   */
  hx_scalar hxscale;
  real fscale, f0;
  hx_array atmp;
  int inv, d, n;

  /* parse the function argument string. */
  if (!fn_scan_args(argstr, fn_argdef_scale, &f0, &fscale, &inv))
    throw("failed to parse add arguments");

  /* check if an inversion was specified. */
  if (inv)
    fscale = 1.0 / fscale;

  /* store a local copy of the dimension index. */
  d = dim;

  /* check the dimension index. setting no dimension implies adding a
   * real value to the array.
   */
  if (d < 0)
    d = 0;

  /* check the dimension index (upper bound). */
  if (d >= D->nd)
    throw("dimension index %d out of bounds [0,%u)", d, D->nd);

  /* FIXME: implement first-point scaling in fn_execute_scale() */

  /* allocate a temporary scalar for the scaling operation. */
  if (!hx_scalar_alloc(&hxscale, D->array.d))
    throw("failed to allocate scalar multiplication operand");

  /* allocate a temporary duplicate array for the scaling operation. */
  if (!hx_array_copy(&atmp, &D->array))
    throw("failed to allocate duplicate array");

  /* set up the scalar value for multiplication. */
  n = 1 << d;
  hxscale.x[n] = fscale;

  /* perform the addition. */
  if (!hx_array_mul_scalar(&atmp, &hxscale, &D->array))
    throw("failed to scale by scalar value %f(%d)", fscale, n);

  /* free the allocated temporary scalar and array. */
  hx_scalar_free(&hxscale);
  hx_array_free(&atmp);

  /* return success. */
  return 1;
}

