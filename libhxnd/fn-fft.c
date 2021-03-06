
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

/* fn_fft(): applies a complex radix-2 fast Fourier transform.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application, or -1.
 * @args: function argument definition array.
 */
int fn_fft (datum *D, const int dim, const fn_arg *args) {
  /* declare variables to hold argument values. */
  int alt, neg, inv;
  real dir;

  /* get the argument values from the argdef array. */
  if (!fn_args_get_all(args, &alt, &neg, &inv))
    throw("failed to get fft arguments");

  /* check the dimension index. */
  if (dim < 0 || dim >= D->nd)
    throw("dimension index %d out of bounds [0,%u)", dim, D->nd);

  /* handle the 'alternate' option. */
  if (alt && !hx_array_alternate_sign(&D->array, D->dims[dim].k))
    throw("failed to apply sign alternation");

  /* handle the 'negate' option. */
  if (neg && !hx_array_negate_basis(&D->array, D->dims[dim].d))
    throw("failed to negate imaginaries");

  /* determine which direction to run the transform. */
  dir = (inv ? HX_FFT_REVERSE : HX_FFT_FORWARD);

  /* run the back-shift operation before inverse transforms. */
  if (dir == HX_FFT_REVERSE) {
    /* perform the back-shift. */
    if (!hx_array_shift(&D->array, D->dims[dim].k, D->array.sz[dim] / 2))
      throw("failed to half-shift before ifft");
  }

  /* run the fourier transform operation. */
  if (!hx_array_fftfn(&D->array, D->dims[dim].d, D->dims[dim].k, dir))
    throw("failed to perform fourier transform");

  /* run the shift operation after forward transforms. */
  if (dir == HX_FFT_FORWARD) {
    /* perform the forward shift. */
    if (!hx_array_shift(&D->array, D->dims[dim].k, D->array.sz[dim] / 2))
      throw("failed to half-shift after fft");
  }

  /* change the fourier transform status flag in the datum dimension.
   */
  D->dims[dim].ft = (inv ? 0 : 1);

  /* return success. */
  return 1;
}

