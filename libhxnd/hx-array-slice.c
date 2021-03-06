
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

/* hx_array_slicer(): slice or store a portion of an array based on lower
 * and upper index boundaries. the operation returns its output in a new
 * array, which need not be allocated prior to the slice.
 * @x: pointer to the input array.
 * @y: pointer to the output array.
 * @lower: the lower index bounds.
 * @upper: the upper index bounds.
 * @dir: either HX_ARRAY_SLICER_SLICE or HX_ARRAY_SLICER_STORE.
 */
int hx_array_slicer (hx_array *x, hx_array *y,
                     hx_index lower,
                     hx_index upper,
                     int dir) {
  /* declare a few required variables:
   * @i: general-purpose loop counter.
   * @xycmp: set high if @y needs allocation.
   * @n: number of coefficients per scalar.
   * @ncpy: number of bytes per scalar.
   * @pidx: input array packed linear index.
   * @pidxy: output array packed linear index.
   * @idx: input array index set.
   * @idxy: output array index set.
   * @sznew: output array sizes.
   */
  int i, xycmp, n, ncpy, pidx, pidxy;
  hx_index idx, idxy, sznew;

  /* store the number of coefficients and bytes per hypercomplex scalar. */
  n = x->n;
  ncpy = n * sizeof(real);

  /* allocate three index arrays for use during iteration. */
  idx = hx_index_copy(x->k, lower);
  idxy = hx_index_alloc(x->k);
  sznew = hx_index_alloc(x->k);

  /* check that the index arrays were allocated successfully. */
  if (!idx || !idxy || !sznew)
    throw("failed to allocate %d indices", x->k);

  /* subtract the lower bound from the upper bound. */
  hx_index_diff(x->k, upper, lower, sznew);

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
      xycmp = (xycmp || (y->sz[i] < sznew[i]));
  }

  /* allocate memory for the sliced outputs, but only if the output
   * array does pass the tests performed above that yield @xycmp.
   */
  if (xycmp && dir == HX_ARRAY_SLICER_SLICE &&
      !hx_array_alloc(y, x->d, x->k, sznew))
    throw("failed to allocate slice destination array");

  /* iterate over the array indices to copy. */
  do {
    /* pack the i/o array indices into linear indices. */
    hx_index_pack(x->k, x->sz, idx, &pidx);
    hx_index_pack(x->k, y->sz, idxy, &pidxy);

    /* copy the coefficient memory. */
    switch (dir) {
      /* slice: x ==> y */
      case HX_ARRAY_SLICER_SLICE:
        memcpy(y->x + n * pidxy, x->x + n * pidx, ncpy);
        break;

      /* store: x <== y */
      case HX_ARRAY_SLICER_STORE:
        memcpy(x->x + n * pidx, y->x + n * pidxy, ncpy);
        break;

      /* other: no-op. */
      default:
        break;
    }

    /* increment the output array index. */
    hx_index_incr(x->k, sznew, idxy);
  } while (hx_index_incr_bounded(x->k, lower, upper, idx));

  /* free the allocated index arrays. */
  hx_index_free(idx);
  hx_index_free(idxy);
  hx_index_free(sznew);

  /* return success. */
  return 1;
}

/* hx_array_vector_slicer(): slice or store a linear section from an array,
 * starting and ending at the extents of a given dimension. this is *much*
 * faster for slicing vectors out of nD arrays than hx_array_slice().
 * @x: pointer to the input array.
 * @y: pointer to the output array.
 * @k: array dimension to slice along.
 * @loc: off-dimension slice origin.
 * @dir: either HX_ARRAY_SLICER_SLICE or HX_ARRAY_SLICER_STORE.
 */
int hx_array_vector_slicer (hx_array *x, hx_array *y,
                            int k, int loc, int dir) {
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
  if ((y->d != x->d || y->k != 1 || y->sz[0] < n) &&
      !hx_array_alloc(y, x->d, 1, &n))
    throw("failed to allocate slice destination array");

  /* copy the scalar values into the (vector) output array. */
  for (i = 0, idx = loc; i < n; i++, idx += stride) {
    /* copy the coefficient memory. */
    switch (dir) {
      /* slice: x ==> y */
      case HX_ARRAY_SLICER_SLICE:
        memcpy(y->x + y->n * i, x->x + x->n * idx, ncpy);
        break;

      /* store: x <== y */
      case HX_ARRAY_SLICER_STORE:
        memcpy(x->x + x->n * idx, y->x + y->n * i, ncpy);
        break;

      /* other: no-op. */
      default:
        break;
    }
  }

  /* return success. */
  return 1;
}

/* hx_array_matrix_slicer(): slice or store a planar section from an array,
 * similar to hx_array_vector_slicer().
 * @x: pointer to the input array.
 * @y: pointer to the output array.
 * @k1: first slice array dimension.
 * @k2: second slice array dimension.
 * @loc: off-dimension slice origin.
 * @dir: either HX_ARRAY_SLICER_SLICE or HX_ARRAY_SLICER_STORE.
 */
int hx_array_matrix_slicer (hx_array *x, hx_array *y,
                            int k1, int k2, int loc,
                            int dir) {
  /* declare a few required variables:
   * @i, @j: loop counters to use during slicing and index setup.
   * @n: sizes of the sliced dimensions, output array.
   * @ncpy: number of bytes per hypercomplex scalar value.
   * @idx: linear index of the input array.
   * @idxy: linear index of the output array.
   * @stride: strides for each slice dimension.
   */
  int i, j, n[2], ncpy, idx, idxy, stride[2];

  /* check that the slice dimensions are in order. */
  if (k1 >= k2)
    throw("slice dimensions (%d,%d) out of order", k1, k2);

  /* check that the slice dimensions are in bounds. */
  if (k1 < 0 || k2 < 0 || k1 >= x->k || k2 >= x->k)
    throw("slice dimensions (%d,%d) out of bounds [0,%d)U[0,%d)",
          k1, k2, x->k, x->k);

  /* compute the number of bytes per scalar
   * and the sizes of the output array.
   */
  ncpy = x->n * sizeof(real);
  n[0] = x->sz[k1];
  n[1] = x->sz[k2];

  /* compute the stride along the first slice dimension. */
  for (i = 0, stride[0] = 1; i < k1; i++)
    stride[0] *= x->sz[i];

  /* compute the stride along the second slice dimension. */
  for (j = 0, stride[1] = 1; j < k2; j++)
    stride[1] *= x->sz[j];

  /* adjust the second stride to account for incrementation in the
   * first slice dimension.
   */
  stride[1] -= stride[0] * (n[0] - 1);

  /* allocate the output array, if its configuration does not match
   * the required configuration.
   */
  if ((y->d != x->d || y->k != 2 || y->sz[0] != n[0] || y->sz[1] != n[1]) &&
      !hx_array_alloc(y, x->d, 2, n))
    throw("failed to allocate slice destination array");

  /* arrays are column-major, so vary the larger index slower. */
  for (j = 0, idx = loc; j < n[1]; j++, idx += stride[1]) {
    /* then vary the smaller index faster inside the slow loop. */
    for (i = 0; i < n[0]; i++, idx += (i < n[0] ? stride[0] : 0)) {
      /* compute the matrix coefficient index. */
      idxy = i + j * n[0];

      /* copy the coefficient memory. */
      switch (dir) {
        /* slice: x ==> y */
        case HX_ARRAY_SLICER_SLICE:
          memcpy(y->x + y->n * idxy, x->x + x->n * idx, ncpy);
          break;

        /* store: x <== y */
        case HX_ARRAY_SLICER_STORE:
          memcpy(x->x + x->n * idx, y->x + y->n * idxy, ncpy);
          break;

        /* other: no-op. */
        default:
          break;
      }
    }
  }

  /* return success. */
  return 1;
}

/* hx_array_sched_slicer(): slice or store a scheduled set of points from
 * an array.
 * @x: pointer to the input array.
 * @y: pointer to the output array.
 * @off: array slice offset index.
 * @n: array slice output point count.
 * @sched: array of off-dimension slice origins.
 * @dir: either HX_ARRAY_SLICER_SLICE or HX_ARRAY_SLICER_STORE.
 */
int hx_array_sched_slicer (hx_array *x, hx_array *y,
                           int off, int n, hx_index sched,
                           int dir) {
  /* declare required variables. */
  int idx, idxy, ncpy;

  /* compute the number of bytes per scalar. */
  ncpy = x->n * sizeof(real);

  /* allocate the output array, if its configuration does not match
   * the required configuration.
   */
  if ((y->d != x->d || y->k != 1 || y->sz[0] != n) &&
      !hx_array_alloc(y, x->d, 1, &n))
    throw("failed to allocate slice destination array");

  /* loop over the points in the schedule. */
  for (idxy = 0; idxy < n; idxy++) {
    /* compute the input matrix coefficient index. */
    idx = off + sched[idxy];

    /* copy the coefficient memory. */
    switch (dir) {
      /* slice: x ==> y */
      case HX_ARRAY_SLICER_SLICE:
        memcpy(y->x + y->n * idxy, x->x + x->n * idx, ncpy);
        break;

      /* store: x <== y */
      case HX_ARRAY_SLICER_STORE:
        memcpy(x->x + x->n * idx, y->x + y->n * idxy, ncpy);
        break;

      /* other: no-op. */
      default:
        break;
    }
  }

  /* return success. */
  return 1;
}

