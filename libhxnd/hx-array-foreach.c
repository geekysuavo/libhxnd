
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

/* include the n-dimensional math header. */
#include <hxnd/hx.h>

/* hx_array_foreach_vector(): perform an operation on each vector of a
 * hypercomplex array using a standardized callback-based scheme.
 * @x: pointer to the array to manipulate.
 * @k: topological dimension (mode) of the vectors.
 * @fn: per-vector callback function pointer.
 */
int hx_array_foreach_vector (hx_array *x, int k,
                             hx_array_foreach_cb fn, ...) {
  /* declare a few required variables:
   * @idx: the index array for the current position in @x.
   * @pidx: the packed linear index for the current position in @x.
   * @y: hypercomplex array holding the currently sliced vector values.
   * @vl: the variable argument list passed to each callback invocation.
   */
  int pidx, szk, slice;
  hx_index idx;
  hx_array y;
  va_list vl;

  /* check that the dimension index is in bounds. */
  if (k < 0 || k >= x->k)
    throw("dimension index %d is out of bounds [0,%d)", k, x->k);

  /* retrieve the size of the dimension under operation. */
  szk = x->sz[k];

  /* allocate an index array. */
  idx = hx_index_alloc(x->k);

  /* check that allocation was successful. */
  if (!idx)
    throw("failed to allocate %d indices", x->k);

  /* allocate a temporary array to store each sliced vector. */
  if (!hx_array_alloc(&y, x->d, 1, &szk))
    throw("failed to allocate slice (%d, 1)-array", x->d);

  /* initialize the linear indices. */
  pidx = slice = 0;

  /* iterate over the elements of the array. */
  do {
    /* pack the index array into a linear index. */
    hx_index_pack(x->k, x->sz, idx, &pidx);

    /* slice the currently indexed vector from the array. */
    if (!hx_array_slice_vector(x, &y, k, pidx))
      throw("failed to slice vector %d", slice);

    /* initialize the variable arguments list. */
    va_start(vl, fn);

    /* execute the callback function. */
    if (!fn(x, &y, idx, pidx, &vl))
      throw("failed to execute callback %d", slice);

    /* free the variable arguments list. */
    va_end(vl);

    /* store the modified sliced vector back into the array. */
    if (!hx_array_store_vector(x, &y, k, pidx))
      throw("failed to store vector %d", slice);

    /* increment the slice index. */
    slice++;
  } while (hx_index_skip(x->k, x->sz, idx, k));

  /* free the temporary array. */
  hx_array_free(&y);

  /* free the multidimensional index. */
  hx_index_free(idx);

  /* return success. */
  return 1;
}

/* hx_array_foreach_matrix(): perform an operation on each plane of a
 * hypercomplex array using a standardized callback-based scheme.
 * @x: pointer to the array to manipulate.
 * @k1: first topological dimension of the planes.
 * @k2: second topological dimension of the planes.
 * @fn: per-matrix callback function pointer.
 */
int hx_array_foreach_matrix (hx_array *x, int k1, int k2,
                             hx_array_foreach_cb fn, ...) {
  /* declare a few required variables:
   * @kl: smaller of the two passed array indices.
   * @ku: larger of the two passed array indices.
   * @idx: the multidimensional for the current position in @x.
   * @pidx: the packed linear index for the current position in @x.
   * @sz: the sizes of each sliced submatrix.
   * @mask: array of masked dimensions to skip during incrementation.
   * @slice: index for counting each sliced matrix.
   * @y: hypercomplex array holding the currently sliced matrix values.
   * @vl: the variable argument list passed to each callback invocation.
   */
  int kl, ku, pidx, slice;
  hx_index idx, sz, mask;
  hx_array y;
  va_list vl;

  /* check that the first dimension index is in bounds. */
  if (k1 < 0 || k1 >= x->k)
    throw("first dimension index %d out of bounds [0,%d)", k1, x->k);

  /* check that the second dimension index is in bounds. */
  if (k2 < 0 || k2 >= x->k)
    throw("second dimension index %d out of bounds [0,%d)", k2, x->k);

  /* allocate temporary index and size arrays. */
  mask = hx_index_alloc(x->k);
  idx = hx_index_alloc(x->k);
  sz = hx_index_alloc(2);

  /* check that allocation was successful. */
  if (!idx || !mask || !sz)
    throw("failed to allocate %d indices", 2 * x->k + 2);

  /* sort the dimension indices. */
  kl = (k1 < k2 ? k1 : k2);
  ku = (k1 < k2 ? k2 : k1);

  /* store the array dimension sizes in the submatrix size array. */
  sz[0] = x->sz[kl];
  sz[1] = x->sz[ku];

  /* allocate a temporary array to store each sliced matrix. */
  if (!hx_array_alloc(&y, x->d, 2, sz))
    throw("failed to allocate slice (%d, 2)-array", x->d);

  /* initialize the linear indices. */
  pidx = slice = 0;

  /* initialize the mask array. */
  mask[kl] = mask[ku] = 1;

  /* iterate over the elements of the array. */
  do {
    /* pack the index array into a linear index. */
    hx_index_pack(x->k, x->sz, idx, &pidx);

    /* slice the currently indexed matrix from the array. */
    if (!hx_array_slice_matrix(x, &y, kl, ku, pidx))
      throw("failed to slice matrix %d", slice);

    /* initialize the variable arguments list. */
    va_start(vl, fn);

    /* execute the callback function. */
    if (!fn(x, &y, idx, pidx, &vl))
      throw("failed to execute callback %d", slice);

    /* free the variable arguments list. */
    va_end(vl);

    /* store the modified sliced matrix back into the array. */
    if (!hx_array_store_matrix(x, &y, kl, ku, pidx))
      throw("failed to store matrix %d", slice);

    /* increment the slice index. */
    slice++;
  } while (hx_index_incr_mask(x->k, x->sz, idx, mask));

  /* free the temporary array. */
  hx_array_free(&y);

  /* free the temporary multidimensional indices. */
  hx_index_free(mask);
  hx_index_free(idx);
  hx_index_free(sz);

  /* return success. */
  return 1;
}

/* hx_array_projector(): compute a hypercomplex scalar value for each vector
 * of a hypercomplex multidimensional array and return the values in a new
 * (uncompacted) array.
 * @x: pointer to the array to project.
 * @k: topological dimension to project through.
 * @fn: per-vector projection function pointer.
 * @xp: pointer to the projected output array.
 */
int hx_array_projector (hx_array *x, int k, hx_array_projector_cb fn,
                        hx_array *xp) {
  /* declare a few required variables:
   * @idx: the multidimensional index for the current position in @x.
   * @pidx: the packed linear index for the current position in @y.
   * @szk: the topological size along the slice dimension.
   * @sznew: topological size of the 'projected' output array.
   * @slice: linear slice index used for counting progress.
   * @y: hypercomplex array holding the currently sliced vector values.
   */
  int pidx, szk, slice;
  hx_index idx, sznew;
  hx_array y;

  /* check that the dimension index is in bounds. */
  if (k < 0 || k >= x->k)
    throw("dimension index %d is out of bounds [0,%d)", k, x->k);

  /* duplicate the source size array, allocate an index array. */
  sznew = hx_index_copy(x->k, x->sz);
  idx = hx_index_alloc(x->k);

  /* check that allocation was successful. */
  if (!sznew || !idx)
    throw("failed to allocate pair of %d indices", x->k);

  /* retrieve the size of the dimension under operation and
   * adjust the output size array.
   */
  szk = sznew[k];
  sznew[k] = 1;

  /* allocate a temporary array to store each sliced vector. */
  if (!hx_array_alloc(&y, x->d, 1, &szk))
    throw("failed to allocate slice (%d, 1)-array", x->d);

  /* check if the output array requires allocation. */
  if (xp->d != x->d || xp->k != x->k ||
      hx_index_cmp(xp->k, xp->sz, sznew)) {
    /* yes. allocate the output array. */
    if (!hx_array_alloc(xp, x->d, x->k, sznew))
      throw("failed to allocate projection array");
  }

  /* initialize the linear indices. */
  pidx = slice = 0;

  /* iterate over the elements of the array. */
  do {
    /* pack the index array into an output linear index. */
    hx_index_pack(x->k, sznew, idx, &pidx);

    /* slice the currently indexed vector from the array. */
    if (!hx_array_slice_vector(x, &y, k, pidx))
      throw("failed to slice vector %d", slice);

    /* execute the callback function. */
    if (!fn(&y, xp->x + pidx * xp->n))
      throw("failed to execute callback %d", slice);

    /* increment the slice index. */
    slice++;
  } while (hx_index_skip(x->k, x->sz, idx, k));

  /* free the temporary array. */
  hx_array_free(&y);

  /* free the allocated multidimensional indices. */
  hx_index_free(sznew);
  hx_index_free(idx);

  /* return success. */
  return 1;
}

