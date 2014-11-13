
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

/* fn_argdef_fft: define all accepted arguments for the 'fft' function.
 */
static const fn_args fn_argdef_fft[] = {
  { "alternate", FN_ARGTYPE_BOOL, "0" },
  { "negate",    FN_ARGTYPE_BOOL, "0" },
  { "inverse",   FN_ARGTYPE_BOOL, "0" },
  { NULL, '\0', NULL }
};

/* fn_execute_fft(): applies a complex radix-2 fast Fourier transform.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application.
 * @args: function argument string.
 */
int fn_execute_fft (datum *D, const int dim, const char *argstr) {
  /* declare variables to hold argument values.
   */
  int alt, neg, inv;
  real dir;

  /* parse the function argument string. */
  if (!fn_scan_args(argstr, fn_argdef_fft, &alt, &neg, &inv))
    throw("failed to parse fft arguments");

  /* check the dimension index. */
  if (dim < 0 || dim >= D->nd)
    throw("dimension index %d out of bounds [0,%u)", dim, D->nd);

  /* FIXME: handle 'alternate' option. */

  /* FIXME: handle 'negate' option. */

  /* determine which direction to run the transform. */
  dir = (inv ? HX_FFT_REVERSE : HX_FFT_FORWARD);

  /* run the fourier transform operation. */
  if (!hx_array_fftfn(&D->array, dim, dim, dir))
    throw("failed to perform fourier transform");

  /* change the fourier transform status flag in the datum dimension.
   */
  D->dims[dim].ft = 1;

  /* return success. */
  return 1;
}

