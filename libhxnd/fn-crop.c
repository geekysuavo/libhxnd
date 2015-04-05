
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

/* fn_crop(): selects a subregion of one dimension of a datum structure.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application, or -1.
 * @args: function argument definition array.
 */
int fn_crop (datum *D, const int dim, const fn_arg *args) {
  /* declare variables to hold argument values:
   * @flo: lower-bound cropping frequency.
   * @fhi: upper-bound cropping frequency.
   * @ppm: whether the crop locations are in ppm or not.
   * @hz: whether the crop locations are in hz or not.
   */
  real flo, fhi;
  int ppm, hz;

  /* declare a few required variables:
   * @ldim: local datum dimension index.
   * @k: topological dimension index of the crop dimension.
   * @szk: topological size of the crop dimension.
   * @car: carrier frequency of the crop dimension, in mhz.
   * @sw: spectral width of the crop dimension, in hertz.
   * @off: offset frequency of the crop dimension, in hertz.
   */
  int ldim, k, szk, szknew, ilo, ihi;
  real car, sw, off, swnew, offnew;
  hx_index sznew;

  /* get the argument values from the argdef array. */
  if (!fn_args_get_all(args, &flo, &fhi, &ppm, &hz))
    throw("failed to get crop arguments");

  /* store the dimensionality into a local variable. */
  ldim = dim;
  if (ldim < 0)
    ldim = 0;

  /* check the dimension index. */
  if (ldim >= D->nd)
    throw("dimension index %d out of bounds [0,%d)", ldim, D->nd);

  /* check that the dimension is frequency-domain. */
  if (!D->dims[ldim].ft)
    throw("dimension %d is not frequency-domain", ldim);

  /* check that no more than one unit was specified. */
  if (ppm && hz)
    throw("multiple unit options set");

  /* get the array dimensionalities. */
  k = D->dims[ldim].k;
  szk = D->array.sz[k];

  /* extract the carrier frequency, spectral width and offset. */
  car = (D->dims[ldim].carrier == 0.0 ? 1.0 : D->dims[ldim].carrier);
  sw = (D->dims[ldim].width == 0.0 ? 1.0 : D->dims[ldim].width);
  off = D->dims[ldim].offset;

  /* check if the crop bound frequencies are in ppm. */
  if (ppm) {
    /* convert the frequencies from ppm to hertz. */
    flo *= car;
    fhi *= car;
  }

  /* check if the crop bound frequencies are now in hertz. */
  if (ppm || hz) {
    /* convert the frequencies from hertz to normalized units. */
    flo = ((flo - off) / sw) + 0.5;
    fhi = ((fhi - off) / sw) + 0.5;
  }

  /* quantize the lower bound frequency by the digital resolution. */
  flo = floor(flo * (real) (szk - 1));
  ilo = (int) flo;
  flo /= (real) (szk - 1);

  /* quantize the upper bound frequency as well. */
  fhi = ceil(fhi * (real) (szk - 1));
  ihi = (int) fhi;
  fhi /= (real) (szk - 1);

  /* compute the new spectral point count along the crop dimension. */
  szknew = ihi - ilo + 1;

  /* compute the new spectral offset and width. */
  offnew = (sw / 2.0) * (flo + fhi - 1.0) + off;
  swnew = sw * fabs(fhi - flo);

  /* allocate a new size array for crop resizing. */
  sznew = hx_index_copy(D->array.k, D->array.sz);

  /* check that the size array was allocated. */
  if (!sznew)
    throw("failed to copy %d indices", D->array.k);

  /* store the new point count along the crop dimension. */
  sznew[k] = szknew;

  /* check that the lower bound frequency is in bounds. */
  if (flo < 0.0 || flo > 1.0)
    throw("lower bound frequency %.3f out of bounds [0,1]", flo);

  /* check that the upper bound frequency is in bounds. */
  if (fhi < 0.0 || fhi > 1.0)
    throw("upper bound frequency %.3f out of bounds [0,1]", fhi);

  /* check that the lower bound is less than the upper bound. */
  if (flo >= fhi)
    throw("upper bound frequency must exceed the lower bound");

  /* shift the array contents along the crop dimension. */
  if (!hx_array_shift(&D->array, k, -ilo))
    throw("failed to perform left shift");

  /* crop the datum array to contain only the desired data points. */
  if (!hx_array_resize(&D->array, D->array.d, D->array.k, sznew))
    throw("failed to resize core array");

  /* store the newly computed spectral width and offset values. */
  D->dims[ldim].offset = offnew;
  D->dims[ldim].width = swnew;
  D->dims[ldim].sz = szknew;

  /* free the size array. */
  hx_index_free(sznew);

  /* return success. */
  return 1;
}

