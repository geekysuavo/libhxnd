
/* hxnd: A framework for n-dimensional hypercomplex calculations for NMR.
 * Copyright (C) 2014-2015  Bradley Worley  <geekysuavo@gmail.com>.
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

/* include the n-dimensional math header. */
#include <hxnd/hx.h>

/* hx_array_ffm1d(): perform a set of one-dimensional fast forward
 * maximum entropy reconstructions over a two-dimensional array.
 * see hx_array_ffm() for details.
 */
int hx_array_ffm1d (hx_array *x, hx_index dx, hx_index kx,
                    int nsched, hx_index sched, int niter,
                    hx_entropy_functional f,
                    hx_entropy_functional df) {
  /* FIXME: implement hx_array_ffm1d() */

  /* return success. */
  return 1;
}

/* hx_array_ffmnd(): perform a set of multidimensional fast forward
 * maximum entropy reconstructions over an array having at least
 * three dimensions.
 * see hx_array_ffm() for details.
 */
int hx_array_ffmnd (hx_array *x, hx_index dx, hx_index kx,
                    int dsched, int nsched, hx_index sched, int niter,
                    hx_entropy_functional f,
                    hx_entropy_functional df) {
  /* FIXME: implement hx_array_ffm2d() */

  /* return success. */
  return 1;
}

/* hx_array_ffm(): use the fast forward maximum entropy algorithm to
 * reconstruct nonuniformly subsampled time-domain data in a hypercomplex
 * array.
 * @x: pointer to the array to reconstruct.
 * @dx: array of algebraic dimension indices in @x.
 * @kx: array of topological dimension indices in @x.
 * @dsched: number of reconstruction dimensions.
 * @nsched: number of sampled @dsched-dimensional points.
 * @niter: maximum number of iterations to perform.
 * @type: entropy functional type to utilize..
 */
int hx_array_ffm (hx_array *x, hx_index dx, hx_index kx,
                  int dsched, int nsched, hx_index sched,
                  int niter, enum hx_entropy_type type) {
  /* declare a few required variables:
   * @f: function pointer for the entropy functional.
   * @df: function pointer for the entropy derivative.
   */
  hx_entropy_functional f, df;

  /* ensure the schedule is allocated and properly sized. */
  if (!sched || dsched < 1 || nsched < 1)
    throw("invalid schedule configuration (%dx%d)", dsched, nsched);

  /* ensure the iteration count is in bounds. */
  if (niter < 1)
    throw("iteration count %d out of bounds [1,inf)", niter);

  /* retrieve the entropy functional function pointers. */
  if (!hx_entropy_get_functionals(type, &f, &df))
    throw("failed to retrieve entropy functionals");

  /* determine which reconstruction function to use. */
  if (dsched == 1) {
    /* execute the one-dimensional function. */
    if (!hx_array_ffm1d(x, dx, kx, nsched, sched, niter, f, df))
      throw("failed to execute one-dimensional reconstruction");
  }
  else if (dsched > 1) {
    /* execute the multidimensional function. */
    if (!hx_array_ffmnd(x, dx, kx, dsched, nsched, sched, niter, f, df))
      throw("failed to execute n-dimensional reconstruction");
  }
  else {
    /* throw an exception. */
    throw("invalid schedule configuration (%dx%d)", dsched, nsched);
  }

  /* return success. */
  return 1;
}

