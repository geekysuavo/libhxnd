
/* ndmath: A framework for n-dimensional hypercomplex calculations for NMR.
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

/* ensure once-only inclusion. */
#ifndef __HXND_HX_ARRAY_H__
#define __HXND_HX_ARRAY_H__

/* hx_array: data type for nD arrays of hypercomplex nD numbers.
 *
 * all the above information still applies to arrays, but an extra layer of
 * complication is added to arrays, namely: they are arrays. instead of just
 * using literal arrays of hx_scalar structures, the array data is squashed
 * together into a single C array of real scalars.
 *
 * as it turns out, this squashing not only aids in memory compactness of
 * multidimensional datasets, but also in communicating array data between
 * the opencl host and device layers.
 */
typedef struct {
  /* d: dimensionality of the hypercomplex space.
   * n: number of coefficients per hypercomplex value (2**d).
   */
  int d, n;

  /* k: dimensionality of the array.
   * len: number of total array coefficients.
   * sz: array sizes along each array dimension.
   */
  int k, len, *sz;

  /* real coefficients. */
  real *x;

  /* multiplication table. this is also a shared table. */
  hx_algebra tbl;
}
hx_array;

/* function declarations: */

int hx_array_alloc (hx_array *x, int d, int k, int *sz);

int hx_array_free (hx_array *x);

int hx_array_set_coeff (hx_array *x, int di, real value, ...);

int hx_array_print (hx_array *x, const char *fname);

int hx_array_save (hx_array *x, const char *fname);

int hx_array_load (hx_array *x, const char *fname);

int hx_array_is_vector (hx_array *x);

int hx_array_is_matrix (hx_array *x);

int hx_array_is_cube (hx_array *x);

int hx_array_resize (hx_array *x, int d, int k, int *sz);

#define hx_array_real(x) \
    hx_array_resize((x), 0, (x)->k, (x)->sz);

int hx_array_reshape (hx_array *x, int k, int *sz);

int hx_array_repack (hx_array *x, int ndiv);

int hx_array_slice (hx_array *x, hx_array *y, int *lower, int *upper);

int hx_array_slice_vector (hx_array *x, hx_array *y, int k, int *loc);

int hx_array_store_vector (hx_array *x, hx_array *y, int k, int *loc);

#endif /* __HXND_HX_ARRAY_H__ */

