
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

/* include the filter function header. */
#include <hxnd/hx-filter.h>

/* fn_filter(): applies a fir filter to a dimension of a datum structure.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application, or -1.
 * @args: function argument definition array.
 */
int fn_filter (datum *D, const int dim, const fn_arg *args) {
  /* declare variables to hold argument values:
   * @order: fir filter order.
   * @ppm: whether the transitions are in ppm or not.
   * @hz: whether the transitions are in hz or not.
   * @flo: low-pass transition frequency.
   * @fhi: high-pass transition frequency.
   */
  int order, ppm, hz;
  real flo, fhi;

  /* declare a few required variables:
   * @sw: spectral width of the filtered dimension, in hertz.
   * @car: carrier frequency of the filtered dimension, in mhz.
   * @f0: normalized center frequency.
   * @ft: normalized transition frequency.
   * @b: array of windowed fir filter coefficients.
   * @ldim: local dimension index for filtering.
   * @inv: whether to invert the filter.
   */
  real sw, car, off, f0, ft;
  int ldim, d, k, szk, inv;
  hx_array b, ph;

  /* get the argument values from the argdef array. */
  if (!fn_args_get_all(args, &order, &flo, &fhi, &ppm, &hz))
    throw("failed to get filter arguments");

  /* store the dimensionality into a local variable. */
  ldim = dim;
  if (ldim < 0)
    ldim = 0;

  /* check the dimension index. */
  if (ldim >= D->nd)
    throw("dimension index %d out of bounds [0,%d)", ldim, D->nd);

  /* check that the dimension is time-domain. */
  if (D->dims[ldim].ft || !D->dims[ldim].cx)
    throw("dimension %d is not complex time-domain", ldim);

  /* check that the filter order is even. */
  if (order < 2 || order % 2)
    throw("filter order must be even");

  /* check that no more than one unit was specified. */
  if (ppm && hz)
    throw("multiple unit options set");

  /* get the array dimensionalities. */
  d = D->dims[ldim].d;
  k = D->dims[ldim].k;
  szk = D->array.sz[k];

  /* extract the carrier frequency, spectral width and offset. */
  car = (D->dims[ldim].carrier == 0.0 ? 1.0 : D->dims[ldim].carrier);
  sw = (D->dims[ldim].width == 0.0 ? 1.0 : D->dims[ldim].width);
  off = D->dims[ldim].offset;

  /* check if the transition frequencies are in ppm. */
  if (ppm) {
    /* convert the transitions from ppm to hertz. */
    flo *= car;
    fhi *= car;
  }

  /* check if the transition frequencies are now in hertz. */
  if (ppm || hz) {
    /* convert the transitions from hertz to normalized frequency. */
    flo = ((flo - off) / sw) + 0.5;
    fhi = ((fhi - off) / sw) + 0.5;
  }

  /* check that the low-pass transition frequency is in bounds. */
  if (isfinite(flo) && (flo < 0.0 || flo > 1.0))
    throw("low-pass frequency %.3f out of bounds [0,1]", flo);

  /* check that the high-pass transition frequency is in bounds. */
  if (isfinite(fhi) && (fhi < 0.0 || fhi > 1.0))
    throw("high-pass frequency %.3f out of bounds [0,1]", fhi);

  /* initialize the inversion flag. */
  inv = 0;

  /* determine which type of filter to construct. */
  if (isfinite(flo) && isfinite(fhi)) {
    /* check that the filter has nonzero bandwidth. */
    if (flo == fhi)
      throw("unsupported filter: zero bandwidth");

    /* determine whether the filter is band-stop (inv) or not. */
    inv = (flo < fhi);
  }
  else if (isfinite(flo)) {
    /* low-pass. set the high-pass frequency to zero. */
    fhi = 0.0;
  }
  else if (isfinite(fhi)) {
    /* high-pass. set the low-pass frequency to one. */
    flo = 1.0;
  }
  else
    throw("filter parameters not specified");

  /* compute the normalized modulation and transition frequencies. */
  f0 = (flo + fhi) / 2.0 - 0.5;
  ft = fabs(fhi - flo) / 2.0;

  /* compute the filter coefficients. */
  hx_array_init(&b);
  if (!hx_filter_fir_alloc(&b, order, ft, inv))
    throw("failed to construct filter array");

  /* allocate a modulation array. */
  if (!hx_array_alloc(&ph, D->array.d, 1, &szk))
    throw("failed to allocate temporary phasor array");

  /* modulate the data to move into the filter center frequency. */
  if (!hx_array_phasor(&ph, d, 0.0, -f0, 0.0) ||
      !hx_array_mul_vector(&D->array, &ph, k, &D->array))
    throw("failed to apply frequency modulation");

  /* apply the computed filter coefficients. */
  if (!hx_filter_fir(&D->array, k, &b))
    throw("failed to perform filtering");

  /* back-shift the data to remove the linear phase error. */
  if (!hx_array_shift(&D->array, k, -order / 2))
    throw("failed to perform linear phase removal");

  /* demodulate the data back to the proper center frequency. */
  if (!hx_array_phasor(&ph, d, 0.0, f0, 0.0) ||
      !hx_array_mul_vector(&D->array, &ph, k, &D->array))
    throw("failed to undo frequency modulation");

  /* free the allocated temporary arrays. */
  hx_array_free(&ph);
  hx_array_free(&b);

  /* return success. */
  return 1;
}

