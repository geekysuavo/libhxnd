
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

/* include the n-dimensional indexing header. */
#include <hxnd/hx-index.h>

/* hx_index_alloc(): allocate an array of multidimensional indices.
 * @k: the array size.
 */
hx_index hx_index_alloc (int k) {
  /* allocate memory for the array, and return it directly. */
  return (hx_index) calloc (k, sizeof(int));
}

/* hx_index_build(): allocate a new index array and fills it with values.
 * @k: the array size.
 * @...: the array values.
 */
hx_index hx_index_build (int k, ...) {
  /* declare a few required variables. */
  hx_index idx;
  va_list vl;
  int i;

  /* allocate an index array. */
  idx = hx_index_alloc(k);

  /* check that the allocation was successful. */
  if (!idx)
    return NULL;

  /* initialize the variable arguments list. */
  va_start(vl, k);

  /* loop over the array values, pulling from the arguments list. */
  for (i = 0; i < k; i++)
    idx[i] = va_arg(vl, int);

  /* free the variable arguments list and return the array. */
  va_end(vl);
  return idx;
}

/* hx_index_free(): free an allocated multidimension index array.
 * @idx: the index array to free.
 */
void hx_index_free (hx_index idx) {
  /* free the index, if its allocated. */
  if (idx) {
    free(idx);
    idx = NULL;
  }
}

/* hx_index_init(): initialize an allocated array of indices.
 * @k: the size of the array.
 * @idx: the array to initialize.
 */
void hx_index_init (int k, hx_index idx) {
  /* zero the array contents. */
  memset(idx, 0, k * sizeof(int));
}

/* hx_index_copy(): duplicate an array of multidimensional indices.
 * @k: the array size.
 * @idx: the input array.
 */
hx_index hx_index_copy (int k, hx_index idx) {
  /* declare a few required variables. */
  hx_index idxnew;

  /* allocate a new array. */
  idxnew = hx_index_alloc(k);
  if (!idxnew)
    return NULL;

  /* copy the elements from the source array. */
  memcpy(idxnew, idx, k * sizeof(int));

  /* return the allocated, duplicated array. */
  return idxnew;
}

/* hx_index_pack(): pack an array of multidimensional indices into a
 * linear index.
 * @k: the size of the array.
 * @sz: the sizes of each array dimension.
 * @idx: the input array of unpacked indices.
 * @pidx: pointer to the output packed linear index.
 */
void hx_index_pack (int k, hx_index sz, hx_index idx, int *pidx) {
  /* define a few required variables. */
  int ki, stride;

  /* loop over the dimensions of the array. */
  for (ki = 0, *pidx = 0, stride = 1; ki < k; ki++) {
    /* add the next index into the linear index. */
    *pidx += idx[ki] * stride;
    stride *= sz[ki];
  }
}

/* hx_index_unpack(): unpack a linear index into an array of separate
 * multidimensional indices.
 * @k: the size of the array.
 * @sz: the sizes of each array dimension.
 * @idx: the output array of unpacked indices.
 * @pidx: the input packed linear index.
 */
void hx_index_unpack (int k, hx_index sz, hx_index idx, int pidx) {
  /* define a few required variables. */
  int ki, redidx;

  /* loop over the dimensions of the array. */
  for (ki = 0, redidx = pidx; ki < k; ki++) {
    /* extract the current index. */
    idx[ki] = redidx % sz[ki];

    /* reduce the index by the current stride. */
    redidx = (redidx - idx[ki]) / sz[ki];
  }
}

/* hx_index_pack_tiled(): pack a multidimensional tile index and a
 * multidimensional point index into a linear index based on tile
 * size and count.
 * @k: the size of the arrays.
 * @ntile: array of tile counts along each dimension.
 * @sztile: array of tile sizes along each dimension.
 * @idx: input array of point indices inside the tile.
 * @idxt: input array of tile indices.
 * @pidx: pointer to the output packed linear index.
 */
void hx_index_pack_tiled (int k, hx_index ntile, hx_index sztile,
                          hx_index idx, hx_index idxt, int *pidx) {
  /* define a few required variables. */
  int ki, idxki, stride;

  /* loop over the dimensions of the array. */
  for (ki = 0, *pidx = 0, stride = 1; ki < k; ki++) {
    /* add the next index into the linear index. */
    idxki = idx[ki] + idxt[ki] * sztile[ki];
    *pidx += idxki * stride;

    /* scale the stride for the next dimension offset. */
    stride *= ntile[ki] * sztile[ki];
  }
}

/* hx_index_incr(): increment a multidimensional index. the function returns
 * '1' except when the increment operation causes all indices to roll over
 * to zero at once (i.e. at the array end).
 *
 * this function is meant to be used as the condition inside
 * a do{}while() loop.
 *
 * @k: the size of the array.
 * @sz: the sizes of each array dimension.
 * @idx: the array of unpacked indices to increment.
 */
int hx_index_incr (int k, hx_index sz, hx_index idx) {
  /* declare a few required variables:
   * @ki: the dimension loop counter.
   * @roundtrip: whether we've made it back to (0,0,0,...)
   */
  int ki, roundtrip;

  /* initialize the round trip indicator. */
  roundtrip = 0;

  /* loop over the dimensions of the array. */
  for (ki = 0; ki < k; ki++) {
    /* increment the current index. */
    idx[ki]++;

    /* check if the current index has overflowed. */
    if (idx[ki] >= sz[ki]) {
      /* reset the index. */
      idx[ki] = 0;
    }
    else {
      /* break the loop. */
      break;
    }

    /* check if a round trip was made. */
    if (ki == k - 1)
      roundtrip = 1;
  }

  /* return whether more indices are available. */
  return !roundtrip;
}

/* hx_index_decr(): decrement a multidimensional index. the function returns
 * '1' except when the decrement operation causes all indices to roll over
 * past zero at once (i.e. at the array start).
 *
 * calling this function with an all-zero index will initialize its contents
 * to the array end. this function is meant to be used as the condition inside
 * a do{}while() loop.
 *
 * @k: the size of the array.
 * @sz: the sizes of each array dimension.
 * @idx: the array of unpacked indices to decrement.
 */
int hx_index_decr (int k, hx_index sz, hx_index idx) {
  /* declare a few required variables:
   * @ki: the dimension loop counter.
   * @allzero: whether we've made it back to (0,0,0,...)
   */
  int ki, allzero;

  /* loop over the dimensions to check if the index is all zeros. */
  for (ki = 0, allzero = 1; ki < k; ki++) {
    /* exit if we encounter a nonzero element. */
    if (idx[ki]) {
      allzero = 0;
      break;
    }
  }

  /* check if the index is all zeros. */
  if (allzero) {
    /* initialize the index to the array end. */
    for (ki = 0; ki < k; ki++)
      idx[ki] = sz[ki] - 1;

    /* return. */
    return 0;
  }

  /* loop over the dimensions of the array. */
  for (ki = 0; ki < k; ki++) {
    /* decrement the current index. */
    idx[ki]--;

    /* check if the current index has underflowed. */
    if (idx[ki] < 0) {
      /* reset the index. */
      idx[ki] = sz[ki] - 1;
    }
    else {
      /* break the loop. */
      break;
    }
  }

  /* return whether more indices are available. */
  return 1;
}

/* hx_index_incr_rev(): increment a multidimensional index, in the opposite
 * sense as hx_index_incr().
 */
int hx_index_incr_rev (int k, hx_index sz, hx_index idx) {
  /* declare a few required variables:
   * @ki: the dimension loop counter.
   * @roundtrip: whether we've made it back to (0,0,0,...)
   */
  int ki, roundtrip;

  /* initialize the round trip indicator. */
  roundtrip = 0;

  /* loop over the dimensions of the array. */
  for (ki = k - 1; ki >= 0; ki--) {
    /* increment the current index. */
    idx[ki]++;

    /* check if the current index has overflowed. */
    if (idx[ki] >= sz[ki]) {
      /* reset the index. */
      idx[ki] = 0;
    }
    else {
      /* break the loop. */
      break;
    }

    /* check if a round trip was made. */
    if (ki == 0)
      roundtrip = 1;
  }

  /* return whether more indices are available. */
  return !roundtrip;
}

/* hx_index_decr_rev(): decrement a multidimensional index, in the opposite
 * sense as hx_index_decr().
 */
int hx_index_decr_rev (int k, hx_index sz, hx_index idx) {
  /* declare a few required variables:
   * @ki: the dimension loop counter.
   * @allzero: whether we've made it back to (0,0,0,...)
   */
  int ki, allzero;

  /* loop over the dimensions to check if the index is all zeros. */
  for (ki = 0, allzero = 1; ki < k; ki++) {
    /* exit if we encounter a nonzero element. */
    if (idx[ki]) {
      allzero = 0;
      break;
    }
  }

  /* check if the index is all zeros. */
  if (allzero) {
    /* initialize the index to the array end. */
    for (ki = 0; ki < k; ki++)
      idx[ki] = sz[ki] - 1;

    /* return. */
    return 0;
  }

  /* loop over the dimensions of the array. */
  for (ki = k - 1; ki >= 0; ki--) {
    /* decrement the current index. */
    idx[ki]--;

    /* check if the current index has underflowed. */
    if (idx[ki] < 0) {
      /* reset the index. */
      idx[ki] = sz[ki] - 1;
    }
    else {
      /* break the loop. */
      break;
    }
  }

  /* return whether more indices are available. */
  return 1;
}

/* hx_index_incr_mask(): increment a multidimensional index in a similar
 * fashion to hx_index_incr(), but skip incrementation over all masked
 * indices (mask[k] == 1).
 * @k: the size of the array.
 * @sz: the sizes of each array dimension.
 * @idx: the array of unpacked indices to increment.
 * @mask: the array of index masks to avoid incrementing.
 */
int hx_index_incr_mask (int k, hx_index sz, hx_index idx, hx_index mask) {
  /* declare a few required variables. */
  int ki, roundtrip;

  /* initialize the round trip indicator. */
  roundtrip = 0;

  /* loop over the dimensions of the array. */
  for (ki = 0; ki < k; ki++) {
    /* skip any masked array indices. */
    if (mask[ki])
      continue;

    /* increment the current index. */
    idx[ki]++;

    /* check if the current index has overflowed. */
    if (idx[ki] >= sz[ki]) {
      /* reset the index. */
      idx[ki] = 0;
    }
    else {
      /* break the loop. */
      break;
    }
  }

  /* check if a round trip was made. */
  if (ki == k)
    roundtrip = 1;

  /* return success. */
  return !roundtrip;
}

/* hx_index_incr_bounded(): increment a multidimensional index in a similar
 * fashion to hx_index_incr(), but skip incrementation over all indices not
 * within a set of specified bounds. initialization must be performed before
 * calling this function.
 * @k: the size of the array.
 * @lower: the lower bound of incrementation.
 * @upper: the upper bound of incrementation.
 * @idx: the array of unpacked indices to increment.
 */
int hx_index_incr_bounded (int k, hx_index lower, hx_index upper,
                           hx_index idx) {
  /* declare a few required variables. */
  int ki, roundtrip;

  /* initialize the round trip indicator. */
  roundtrip = 0;

  /* loop over the dimensions of the array. */
  for (ki = 0; ki < k; ki++) {
    /* skip any indices not included in the bounds. */
    if (lower[ki] == upper[ki])
      continue;

    /* increment the current index. */
    idx[ki]++;

    /* check if the current index has overflowed. */
    if (idx[ki] > upper[ki]) {
      /* reset the index. */
      idx[ki] = lower[ki];
    }
    else {
      /* break the loop. */
      break;
    }
  }

  /* check if a round trip was made. */
  if (ki == k)
    roundtrip = 1;

  /* return success. */
  return !roundtrip;
}

/* hx_index_skip(): increment a multidimensional index in a similar fashion
 * to hx_index_incr(), but skip the incrementation over a specified index.
 * @k: the size of the array.
 * @sz: the sizes of each array dimension.
 * @idx: the array of unpacked indices to increment.
 * @kskip: the array index to avoid incrementing.
 */
int hx_index_skip (int k, hx_index sz, hx_index idx, int kskip) {
  /* declare a few required variables. */
  int ki, roundtrip;

  /* initialize the round trip indicator. */
  roundtrip = 0;

  /* loop over the dimensions of the array. */
  for (ki = 0; ki < k; ki++) {
    /* skip the specified array index. */
    if (ki == kskip)
      continue;

    /* increment the current index. */
    idx[ki]++;

    /* check if the current index has overflowed. */
    if (idx[ki] >= sz[ki]) {
      /* reset the index. */
      idx[ki] = 0;
    }
    else {
      /* break the loop. */
      break;
    }
  }

  /* check if a round trip was made. */
  if (ki == k)
    roundtrip = 1;

  /* return success. */
  return !roundtrip;
}

/* hx_index_jump_init(): compute the 'small' and 'large' index strides to
 * use when computing packed linear indices where a single specified
 * index is skipped.
 * @k: the size of the array.
 * @sz: the sizes of each array dimension.
 * @kskip: the array index to avoid incrementing.
 * @ja: pointer to the small stride value.
 * @jb: pointer to the large stride value.
 */
int hx_index_jump_init (int k, hx_index sz, int kskip,
                        int *ja, int *jb, int *jmax) {
  /* declare a required variable. */
  int i;

  /* loop to compute the jump values. */
  for (i = 0, *ja = 1, *jb = sz[0]; i < kskip; i++) {
    *ja *= sz[i];
    *jb *= sz[i + 1];
  }

  /* loop to compute the maximum loop control variable. */
  for (i = 0, *jmax = 1; i < k; i++)
    *jmax *= (i == kskip ? 1 : sz[i]);

  /* return success. */
  return 1;
}

/* hx_index_jump(): compute a linear index from small and large index
 * strides computed by hx_index_jump_init().
 * @j: the for loop control variable, ranges from 0 .. jmax-1.
 * @ja: the small stride value.
 * @jb: the large stride value.
 */
inline int hx_index_jump (int j, int ja, int jb) {
  /* return the computed value. */
  return (jb * (j / ja) + j % ja);
}

/* hx_index_diff(): compute the difference of two multidimensional indices.
 * @k: the number of index elements.
 * @a: the first operand.
 * @b: the second operand.
 * @c: the output index array.
 */
int hx_index_diff (int k, hx_index a, hx_index b, hx_index c) {
  /* declare a required variable. */
  int ki;

  /* compute the element-wise differences. */
  for (ki = 0; ki < k; ki++)
    c[ki] = a[ki] - b[ki];

  /* return success. */
  return 1;
}

/* hx_index_cmp(): compare the values in two multidimensional indices
 * of the same size.
 * @k: the lengths of the two index arrays.
 * @a: the first index array.
 * @b: the second index array.
 */
int hx_index_cmp (int k, hx_index a, hx_index b) {
  /* declare a required variable. */
  int ki;

  /* return if either of the index arrays is unallocated. */
  if (!a || !b)
    return 1;

  /* loop over the indices. */
  for (ki = 0; ki < k; ki++) {
    /* check if the current dimension is a mismatch. */
    if (a[ki] != b[ki])
      return 1;
  }

  /* return identity. */
  return 0;
}

/* hx_index_bounded(): determine whether a multidimensional index is bounded
 * between a (optional) lower index bound and an upper index bound.
 * @k: the index array element count.
 * @idx: the index array in question.
 * @lower: the lower bound array, or NULL for all zeros.
 * @upper: the upper bound array.
 */
int hx_index_bounded (int k, hx_index idx, hx_index lower, hx_index upper) {
  /* declare a few required variables:
   * @ki: index loop counter.
   */
  int ki;

  /* loop over the array elements. */
  for (ki = 0; ki < k; ki++) {
    /* check the lower bound. */
    if ((lower && idx[ki] < lower[ki]) || idx[ki] < 0)
      return 0;

    /* check the upper bound. */
    if (idx[ki] > upper[ki])
      return 0;
  }

  /* the index is within the bounds. */
  return 1;
}

/* hx_index_sort(): sort the values in an index, returning the resulting
 * zero-based indices in their new swapped positions.
 */
int hx_index_sort (int k, hx_index idx) {
  /* declare required variables:
   */
  unsigned int i, j;
  hx_index ord;
  int swp;

  /* allocate a temporary array of indices. */
  ord = hx_index_alloc(k);

  /* ensure the array was allocated. */
  if (!ord)
    return 0;

  /* copy the ordering into the temporary array. */
  memcpy(ord, idx, k * sizeof(int));

  /* initialize the variables in the output array. */
  for (i = 0; i < k; i++)
    idx[i] = i;

  /* loop over the array of dimensions. */
  for (i = 1; i < k; i++) {
    /* set the initial inner loop index. */
    j = i;

    /* loop over the unsorted dimensions. */
    while (j > 0 && ord[j - 1] > ord[j]) {
      /* swap the (j-1) and (j) output indices. */
      swp = idx[j];
      idx[j] = idx[j - 1];
      idx[j - 1] = swp;

      /* swap the ordering of the temporary indices. */
      swp = ord[j];
      ord[j] = ord[j - 1];
      ord[j - 1] = swp;

      /* decrement the inner loop counter. */
      j--;
    }
  }

  /* free the ordering array. */
  free(ord);

  /* return success. */
  return 1;
}

/* hx_index_scheduled(): return a list of packed linear indices locating the
 * sampled elements of an array based on a schedule.
 */
hx_index hx_index_scheduled (int k, hx_index sz, int dsched, int nsched,
                             hx_index sched) {
  /* declare a few required variables. */
  hx_index idx, arr;
  int i, j, nidx;

  /* get the number of scheduled elements. */
  nidx = nsched;

  /* allocate the array of indices and a temporary array. */
  idx = hx_index_alloc(nidx);
  arr = hx_index_alloc(k);

  /* check that allocation was successful. */
  if (!idx || !arr)
    return NULL;

  /* build the array of sampled indices. */
  for (i = 0; i < nidx; i++) {
    /* build the unpacked index array. */
    for (j = 0; j < k; j++)
      arr[j] = sched[i * dsched + j];

    /* pack the index array into a linear index. */
    hx_index_pack(k, sz, arr, idx + i);
  }

  /* free the temporary array and return the indices. */
  hx_index_free(arr);
  return idx;
}

/* hx_index_unscheduled(): return a list of packed linear indices locating
 * the non-sampled elements of an array based on a schedule.
 */
hx_index hx_index_unscheduled (int k, hx_index sz, int dsched, int nsched,
                               hx_index sched) {
  /* declare a few required variables. */
  int i, iadj, j, ntotal, nidx;
  hx_index idx, idxinv;

  /* get the number of total elements. */
  for (i = 0, ntotal = 1; i < k; i++)
    ntotal *= sz[i];

  /* get the number of unscheduled elements. */
  nidx = ntotal - nsched;

  /* allocate the array of indices and a temporary array. */
  idx = hx_index_alloc(nidx);

  /* build the array of scheduled indices. */
  idxinv = hx_index_scheduled(k, sz, dsched, nsched, sched);

  /* check that allocation was successful. */
  if (!idx || !idxinv)
    return NULL;

  /* build the array of unsampled indices. */
  for (i = 0, iadj = 0; i < ntotal; i++) {
    /* search for the index in the sampled array. */
    for (j = 0; j < nsched; j++) {
      /* break if we find the index. */
      if (idxinv[j] == i)
        break;
    }

    /* check if the index has not been sampled. */
    if (j >= nsched) {
      /* store the index and increment the counter. */
      idx[iadj++] = i;
    }
  }

  /* free the scheduled array and return the indices. */
  hx_index_free(idxinv);
  return idx;
}

/* hx_index_printfn(): core function used by hx_index_print() to write the
 * contents of a multidimensional index to standard error.
 * @k: the array size.
 * @idx: the array of indices.
 * @s: the array variable name.
 */
void hx_index_printfn (int k, hx_index idx, const char *s) {
  /* declare a required variable. */
  int ki;

  /* print the variable name. */
  fprintf(stderr, "%s[%d] = (", s, k);

  /* loop over the array elements. */
  for (ki = 0; ki < k; ki++) {
    /* print the index value. */
    fprintf(stderr, "%d", idx[ki]);

    /* print commas to separate values. */
    if (ki < k - 1)
      fprintf(stderr, ", ");
  }

  /* end the line. */
  fprintf(stderr, ")\n");
  fflush(stderr);
}

