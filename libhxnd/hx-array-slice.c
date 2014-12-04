
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

/* hx_array_slicer(): slices a portion of an array based on @lower and @upper
 * index boundaries. the operation returns its output in a new array @y,
 * which should not be allocated prior to the slice.
 * @x: pointer to the input array.
 * @y: pointer to the output array.
 * @lower: the lower index bounds.
 * @upper: the upper index bounds.
 * @dir: either HX_ARRAY_SLICER_SLICE or HX_ARRAY_SLICER_STORE.
 */
int hx_array_slicer (hx_array *x, hx_array *y,
                     int *lower, int *upper,
                     int dir) {
  /* declare a few required variables:
   * @i: general-purpose loop counter.
   * @xycmp: set high if @y needs allocation.
   * @n: number of coefficients per scalar.
   * @ncpy: number of bytes per scalar.
   * @idxi: input array linear index.
   * @idxo: output array linear index.
   * @arri: input array index set.
   * @arro: output array index set.
   * @sznew: output array sizes.
   */
  int i, xycmp, n, ncpy, idxi, idxo, *arri, *arro, *sznew;

  /* store the number of coefficients and bytes per hypercomplex scalar. */
  n = x->n;
  ncpy = n * sizeof(real);

  /* allocate three index arrays for use during iteration. */
  arri = hx_array_index_alloc(x->k);
  arro = hx_array_index_alloc(x->k);
  sznew = hx_array_index_alloc(x->k);

  /* check that the index arrays were allocated successfully. */
  if (!arri || !arro || !sznew)
    throw("failed to allocate %d indices", x->k);

  /* subtract the lower bound from the upper bound. */
  hx_array_index_diff(x->k, upper, lower, sznew);

  /* increment each element of the difference array, resulting in
   * the array of sizes of the sliced portion of the array.
   */
  for (i = 0; i < x->k; i++)
    sznew[i]++;

  /* check that the output array has the same dimensionalities as
   * the input array.
   */
  xycmp = ((x->d != y->d) || (x->k != y->k));

  /* if the above comparisons both evaluated to false, then the following
   * block is safe to execute.
   */
  if (!xycmp) {
    /* check that the output array has the desired size. */
    for (i = 0; i < x->k; i++)
      xycmp = (xycmp || (y->sz[i] != sznew[i]));
  }

  /* allocate memory for the sliced outputs, but only if the output
   * array does pass the tests performed above that yield @xycmp.
   */
  if (xycmp && !hx_array_alloc(y, x->d, x->k, sznew))
    throw("failed to allocate slice destination array");

  /* iterate over the larger (input) array. */
  idxi = 0;
  do {
    /* check if the current input array index is in the slice bounds. */
    if (hx_array_index_bounded(x->k, arri, lower, upper)) {
      /* yes. compute the output array indices. */
      for (i = 0; i < x->k; i++)
        arro[i] = arri[i] - lower[i];

      /* linearize the output indices. */
      hx_array_index_pack(x->k, sznew, arro, &idxo);

      /* copy the coefficient memory. */
      switch (dir) {
        /* slice: x ==> y */
        case HX_ARRAY_SLICER_SLICE:
          memcpy(y->x + n * idxo, x->x + n * idxi, ncpy);
          break;

        /* store: x <== y */
        case HX_ARRAY_SLICER_STORE:
          memcpy(x->x + n * idxi, y->x + n * idxo, ncpy);
          break;

        /* other: no-op. */
        default:
          break;
      }
    }

    /* incremenet the input array linear index. */
    idxi++;
  } while (hx_array_index_incr(x->k, x->sz, arri));

  /* free the allocated index arrays. */
  free(sznew);
  free(arri);
  free(arro);

  /* return success. */
  return 1;
}

/* hx_array_slice_vector(): slice a linear section from an array, starting
 * and ending at the extents of a given dimension @k. this is *much* faster
 * for slicing vectors out of nD arrays than hx_array_slice().
 * @x: pointer to the input array.
 * @y: pointer to the output array.
 * @k: array dimension to slice along.
 * @loc: off-dimension slice origin.
 */
int hx_array_slice_vector (hx_array *x, hx_array *y, int k, int loc) {
  /* declare a few required variables:
   * @i: a loop counter used during slicing and index setup.
   * @n: size of the sliced dimension, length of the output array.
   * @ncpy: number of bytes per hypercomplex scalar value.
   * @idx: linear index of the input array.
   * @stride: stride of the input array.
   */
  int i, n, ncpy, idx, stride;

  /* check that the slice dimension in within bounds. */
  if (k < 0 || k >= x->k)
    throw("slice dimension %d out of bounds [0,%d)", k, x->k);

  /* compute the number of bytes per scalar
   * and the size of the output array.
   */
  ncpy = x->n * sizeof(real);
  n = x->sz[k];

  /* compute the stride for passing along the slice dimension. */
  for (i = 0, stride = 1; i < k; i++)
    stride *= x->sz[i];

  /* allocate the output array, if its configuration does not match
   * the required configuration.
   */
  if ((y->d != x->d || y->k != 1 || y->sz[0] != n) &&
      !hx_array_alloc(y, x->d, 1, &n))
    throw("failed to allocate slice destination array");

  /* copy the scalar values into the (vector) output array. */
  for (i = 0, idx = loc; i < n; i++, idx += stride)
    memcpy(y->x + y->n * i, x->x + x->n * idx, ncpy);

  /* return success. */
  return 1;
}

/* hx_array_store_vector(): store a linear section from an array, essentially
 * the reverse operation of hx_array_slice_vector().
 */
int hx_array_store_vector (hx_array *x, hx_array *y, int k, int loc) {
  /* declare a few required variables. */
  int i, n, ncpy, idx, stride;

  /* check that the slice dimension in within bounds. */
  if (k < 0 || k >= x->k)
    throw("slice dimension %d out of bounds [0,%d)", k, x->k);

  /* compute the number of bytes per scalar
   * and the size of the output array.
   */
  ncpy = x->n * sizeof(real);
  n = x->sz[k];

  /* compute the stride for passing along the slice dimension. */
  for (i = 0, stride = 1; i < k; i++)
    stride *= x->sz[i];

  /* check that the array configurations are correct. */
  if (y->d != x->d || y->k != 1 || y->sz[0] != n)
    throw("source-destination array configuration mismatch");

  /* copy the scalar values into the (vector) output array. */
  for (i = 0, idx = loc; i < n; i++, idx += stride)
    memcpy(x->x + x->n * idx, y->x + y->n * i, ncpy);

  /* return success. */
  return 1;
}

