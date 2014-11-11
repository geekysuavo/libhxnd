
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

/* hx_array_index_alloc(): allocates an array of multidimensional indices.
 */
int *hx_array_index_alloc (int k) {
  /* allocate memory for the array, and return it directly. */
  return (int*) calloc (k, sizeof(int));
}

/* hx_array_index_build(): builds an index array and fills it with values.
 * @k: the array size.
 * @...: the array values.
 */
int *hx_array_index_build (int k, ...) {
  /* declare a few required variables. */
  int i, *arr;
  va_list vl;

  /* allocate an index array. */
  arr = hx_array_index_alloc(k);

  /* check that the allocation was successful. */
  if (!arr) {
    /* raise an error and return nothing. */
    raise("failed to allocate %d indices", k);
    return NULL;
  }

  /* initialize the variable arguments list. */
  va_start(vl, k);

  /* loop over the array values, pulling from the arguments list. */
  for (i = 0; i < k; i++)
    arr[i] = va_arg(vl, int);

  /* free the variable arguments list and return the array. */
  va_end(vl);
  return arr;
}

/* hx_array_index_init(): initializes an allocate array of indices.
 * @arr: the array to initialize.
 * @k: the size of the array.
 */
int hx_array_index_init (int *arr, int k) {
  /* zero the array contents and return success. */
  memset(arr, 0, k * sizeof(int));
  return 1;
}

/* hx_array_index_pack(): packs an array of multidimensional indices into a
 * linear index.
 */
int hx_array_index_pack (int k, int *sz, int *arr, int *pidx) {
  /* define a few required variables. */
  int ki, stride;

  /* loop over the dimensions of the array. */
  for (ki = 0, *pidx = 0, stride = 1; ki < k; ki++) {
    /* add the next index into the linear index. */
    *pidx += arr[ki] * stride;
    stride *= sz[ki];
  }

  /* return success. */
  return 1;
}

/* hx_array_index_unpack(): unpacks a linear index into an array of separate
 * multidimensional indices.
 */
int hx_array_index_unpack (int k, int *sz, int **parr, int idx) {
  /* define a few required variables. */
  int ki, *arr, redidx;

  /* dereference the multidimensional index. */
  arr = *parr;

  /* loop over the dimensions of the array. */
  for (ki = 0, redidx = idx; ki < k; ki++) {
    /* extract the current index. */
    arr[ki] = redidx % sz[ki];

    /* reduce the index by the current stride. */
    redidx = (redidx - arr[ki]) / sz[ki];
  }

  /* return success. */
  return 1;
}

/* hx_array_index_inc(): increments a multidimensional index. the function
 * returns '1' except when the increment operation causes all indices to
 * roll over to zero at once (i.e. at the array end).
 */
int hx_array_index_inc (int k, int *sz, int **parr) {
  /* declare a few required variables:
   * @ki: the dimension loop counter.
   * @arr: the index array.
   * @roundtrip: whether we've made it back to (0,0,0,...)
   */
  int ki, *arr, roundtrip;

  /* initialize the round trip indicator. */
  roundtrip = 0;

  /* dereference the multidimensional index. */
  arr = *parr;

  /* loop over the dimensions of the array. */
  for (ki = 0; ki < k; ki++) {
    /* increment the current index. */
    arr[ki]++;

    /* check if the current index has overflowed. */
    if (arr[ki] >= sz[ki]) {
      /* reset the index. */
      arr[ki] = 0;
    }
    else {
      /* break the loop. */
      break;
    }

    /* check if a round trip was made. */
    if (ki == k - 1)
      roundtrip = 1;
  }

  /* return success. */
  return !roundtrip;
}

/* hx_array_index_diff(): compute the difference of two indices along
 * each dimension.
 * @k: the number of index elements.
 * @a: the first operand.
 * @b: the second operand.
 * @pc: pointer the output array.
 */
int hx_array_index_diff (int k, int *a, int *b, int **pc) {
  /* declare a few required variables:
   * @c: the output index array.
   * @i: a loop counter.
   */
  int *c;
  int i;

  /* dereference the output. */
  c = *pc;

  /* compute the element-wise differences. */
  for (i = 0; i < k; i++)
    c[i] = a[i] - b[i];

  /* return success. */
  return 1;
}

/* hx_array_index_bounded(): determine whether an index @arr is bounded
 * between a (optional) @lower bound and an @upper bound.
 * @k: the index array element count.
 * @arr: the index array in question.
 * @lower: the lower bound array, or NULL for all zeros.
 * @upper: the upper bound array.
 */
int hx_array_index_bounded (int k, int *arr, int *lower, int *upper) {
  /* declare a few required variables:
   * @i: a loop counter.
   */
  int i;

  /* loop over the array elements. */
  for (i = 0; i < k; i++) {
    /* check the lower bound. */
    if ((lower && arr[i] < lower[i]) || arr[i] < 0)
      return 0;

    /* check the upper bound. */
    if (arr[i] > upper[i])
      return 0;
  }

  /* the index is within the bounds. */
  return 1;
}

/* hx_array_index_sort(): sorts the values in @order, returning in @order
 * the resulting zero-based indices in their new swapped positions.
 */
int hx_array_index_sort (int k, int *order) {
  /* declare required variables:
   */
  unsigned int i, j;
  int *zord, swp;

  /* allocate a temporary array of indices. */
  zord = (int*) malloc(k * sizeof(int));
  if (!zord)
    return 0;

  /* copy the ordering into the temporary array. */
  memcpy(zord, order, k * sizeof(int));

  /* initialize the variables in the output array. */
  for (i = 0; i < k; i++)
    order[i] = i;

  /* loop over the array of dimensions. */
  for (i = 1; i < k; i++) {
    /* set the initial inner loop index. */
    j = i;

    /* loop over the unsorted dimensions. */
    while (j > 0 && zord[j - 1] > zord[j]) {
      /* swap the (j-1) and (j) output indices. */
      swp = order[j];
      order[j] = order[j - 1];
      order[j - 1] = swp;

      /* swap the ordering of the temporary indices. */
      swp = zord[j];
      zord[j] = zord[j - 1];
      zord[j - 1] = swp;

      /* decrement the inner loop counter. */
      j--;
    }
  }

  /* free the ordering array. */
  free(zord);

  /* return success. */
  return 1;
}

