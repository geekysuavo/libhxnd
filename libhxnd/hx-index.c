
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

/* hx_array_index_copy(): duplicate an array of multidimensional indices.
 * @k: the array size.
 * @sz: the input array.
 */
int *hx_array_index_copy (int k, int *sz) {
  /* declare a few required variables. */
  int *sznew;

  /* allocate a new array. */
  sznew = hx_array_index_alloc(k);
  if (!sznew)
    return NULL;

  /* copy the elements from the source array. */
  memcpy(sznew, sz, k * sizeof(int));

  /* return the allocated, duplicated array. */
  return sznew;
}

/* hx_array_index_pack(): packs an array of multidimensional indices into a
 * linear index.
 * @k: the size of the array.
 * @sz: the sizes of each array dimension.
 * @arr: the input array of unpacked indices.
 * @pidx: pointer to the output packed linear index.
 */
void hx_array_index_pack (int k, int *sz, int *arr, int *pidx) {
  /* define a few required variables. */
  int ki, stride;

  /* loop over the dimensions of the array. */
  for (ki = 0, *pidx = 0, stride = 1; ki < k; ki++) {
    /* add the next index into the linear index. */
    *pidx += arr[ki] * stride;
    stride *= sz[ki];
  }
}

/* hx_array_index_unpack(): unpacks a linear index into an array of separate
 * multidimensional indices.
 * @k: the size of the array.
 * @sz: the sizes of each array dimension.
 * @arr: the output array of unpacked indices.
 * @idx: the input packed linear index.
 */
void hx_array_index_unpack (int k, int *sz, int *arr, int idx) {
  /* define a few required variables. */
  int ki, redidx;

  /* loop over the dimensions of the array. */
  for (ki = 0, redidx = idx; ki < k; ki++) {
    /* extract the current index. */
    arr[ki] = redidx % sz[ki];

    /* reduce the index by the current stride. */
    redidx = (redidx - arr[ki]) / sz[ki];
  }
}

/* hx_array_index_pack_tiled(): packs a multidimensional tile index and a
 * multidimensional point index into a linear index based on tile size
 * and count.
 * @k: the size of the arrays.
 * @ntile: array of tile counts along each dimension.
 * @sztile: array of tile sizes along each dimension.
 * @arr: input array of point indices inside the tile.
 * @arrt: input array of tile indices.
 * @pidx: pointer to the output packed linear index.
 */
void hx_array_index_pack_tiled (int k, int *ntile, int *sztile,
                                int *arr, int *arrt, int *pidx) {
  /* define a few required variables. */
  int ki, arrki, stride;

  /* loop over the dimensions of the array. */
  for (ki = 0, *pidx = 0, stride = 1; ki < k; ki++) {
    /* add the next index into the linear index. */
    arrki = arr[ki] + arrt[ki] * sztile[ki];
    *pidx += arrki * stride;

    /* scale the stride for the next dimension offset. */
    stride *= ntile[ki] * sztile[ki];
  }
}

/* hx_array_index_incr(): increments a multidimensional index. the function
 * returns '1' except when the increment operation causes all indices to
 * roll over to zero at once (i.e. at the array end).
 *
 * this function is meant to be used as the condition inside
 * a do{}while() loop.
 *
 * @k: the size of the array.
 * @sz: the sizes of each array dimension.
 * @arr: the array of unpacked indices to increment.
 */
int hx_array_index_incr (int k, int *sz, int *arr) {
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

  /* return whether more indices are available. */
  return !roundtrip;
}

/* hx_array_index_decr(): decrements a multidimensional index. the function
 * returns '1' except when the decrement operation causes all indices to
 * roll over past zero at once (i.e. at the array start).
 *
 * calling this function with an all-zero index will initialize its contents
 * to the array end. this function is meant to be used as the condition inside
 * a do{}while() loop.
 *
 * @k: the size of the array.
 * @sz: the sizes of each array dimension.
 * @arr: the array of unpacked indices to decrement.
 */
int hx_array_index_decr (int k, int *sz, int *arr) {
  /* declare a few required variables:
   * @ki: the dimension loop counter.
   * @allzero: whether we've made it back to (0,0,0,...)
   */
  int ki, allzero;

  /* loop over the dimensions to check if the index is all zeros. */
  for (ki = 0, allzero = 1; ki < k; ki++) {
    /* exit if we encounter a nonzero element. */
    if (arr[ki]) {
      allzero = 0;
      break;
    }
  }

  /* check if the index is all zeros. */
  if (allzero) {
    /* initialize the index to the array end. */
    for (ki = 0; ki < k; ki++)
      arr[ki] = sz[ki] - 1;

    /* return. */
    return 0;
  }

  /* loop over the dimensions of the array. */
  for (ki = 0; ki < k; ki++) {
    /* decrement the current index. */
    arr[ki]--;

    /* check if the current index has underflowed. */
    if (arr[ki] < 0) {
      /* reset the index. */
      arr[ki] = sz[ki] - 1;
    }
    else {
      /* break the loop. */
      break;
    }
  }

  /* return whether more indices are available. */
  return 1;
}

/* hx_array_index_incr_rev(): increments a multidimensional index, in the
 * opposite sense as hx_array_index_incr().
 */
int hx_array_index_incr_rev (int k, int *sz, int *arr) {
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
    if (ki == 0)
      roundtrip = 1;
  }

  /* return whether more indices are available. */
  return !roundtrip;
}

/* hx_array_index_decr_rev(): decrements a multidimensional index, in the
 * opposite sense as hx_array_index_decr().
 */
int hx_array_index_decr_rev (int k, int *sz, int *arr) {
  /* declare a few required variables:
   * @ki: the dimension loop counter.
   * @allzero: whether we've made it back to (0,0,0,...)
   */
  int ki, allzero;

  /* loop over the dimensions to check if the index is all zeros. */
  for (ki = 0, allzero = 1; ki < k; ki++) {
    /* exit if we encounter a nonzero element. */
    if (arr[ki]) {
      allzero = 0;
      break;
    }
  }

  /* check if the index is all zeros. */
  if (allzero) {
    /* initialize the index to the array end. */
    for (ki = 0; ki < k; ki++)
      arr[ki] = sz[ki] - 1;

    /* return. */
    return 0;
  }

  /* loop over the dimensions of the array. */
  for (ki = k - 1; ki >= 0; ki--) {
    /* decrement the current index. */
    arr[ki]--;

    /* check if the current index has underflowed. */
    if (arr[ki] < 0) {
      /* reset the index. */
      arr[ki] = sz[ki] - 1;
    }
    else {
      /* break the loop. */
      break;
    }
  }

  /* return whether more indices are available. */
  return 1;
}

/* hx_array_index_skip(): increments a multidimensional index in a similar
 * fashion to hx_array_index_incr(), but skips the incrementation over a
 * specified index.
 * @k: the size of the array.
 * @sz: the sizes of each array dimension.
 * @arr: the array of unpacked indices to increment.
 * @kskip: the array index to avoid incrementing.
 */
int hx_array_index_skip (int k, int *sz, int *arr, int kskip) {
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
  }

  /* check if a round trip was made. */
  if (ki == k)
    roundtrip = 1;

  /* return success. */
  return !roundtrip;
}

/* hx_array_index_jump_init(): computes the 'small' and 'large' index strides
 * to use when computing packed linear indices where a single specified
 * index is skipped.
 * @k: the size of the array.
 * @sz: the sizes of each array dimension.
 * @kskip: the array index to avoid incrementing.
 * @ja: pointer to the small stride value.
 * @jb: pointer to the large stride value.
 */
int hx_array_index_jump_init (int k, int *sz, int kskip,
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

/* hx_array_index_jump(): computes a linear index from small and large index
 * strides computed by hx_array_index_jump_init().
 * @j: the for loop control variable, ranges from 0 .. jmax-1.
 * @ja: the small stride value.
 * @jb: the large stride value.
 */
inline int hx_array_index_jump (int j, int ja, int jb) {
  /* return the computed value. */
  return (jb * (j / ja) + j % ja);
}

/* hx_array_index_diff(): compute the difference of two indices along
 * each dimension.
 * @k: the number of index elements.
 * @a: the first operand.
 * @b: the second operand.
 * @c: the output array.
 */
int hx_array_index_diff (int k, int *a, int *b, int *c) {
  /* declare a required variable. */
  int i;

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

/* hx_array_index_scheduled(): returns a list of packed linear indices
 * locating the sampled elements of an array based on a schedule.
 */
int *hx_array_index_scheduled (int k, int *sz, int dsched, int nsched,
                               int *sched) {
  /* declare a few required variables. */
  int i, j, nidx, *idx, *arr;

  /* get the number of scheduled elements. */
  nidx = nsched;

  /* allocate the array of indices and a temporary array. */
  idx = hx_array_index_alloc(nidx);
  arr = hx_array_index_alloc(k);

  /* check that allocation was successful. */
  if (!idx || !arr) {
    /* raise an exception. */
    raise("failed to allocate %d+%d indices", nidx, k);
    return NULL;
  }

  /* build the array of sampled indices. */
  for (i = 0; i < nidx; i++) {
    /* build the unpacked index array. */
    for (j = 0; j < k; j++)
      arr[j] = sched[i * dsched + j];

    /* pack the index array into a linear index. */
    hx_array_index_pack(k, sz, arr, idx + i);
  }

  /* return the indices. */
  return idx;
}

/* hx_array_index_unscheduled(): returns a list of packed linear indices
 * located the non-sampled elements of an array based on a schedule.
 */
int *hx_array_index_unscheduled (int k, int *sz, int dsched, int nsched,
                                 int *sched) {
  /* declare a few required variables. */
  int i, iadj, j, ntotal, nidx, *idx, *idxinv;

  /* get the number of total elements. */
  for (i = 0, ntotal = 1; i < k; i++)
    ntotal *= sz[i];

  /* get the number of unscheduled elements. */
  nidx = ntotal - nsched;

  /* allocate the array of indices and a temporary array. */
  idx = hx_array_index_alloc(nidx);

  /* build the array of scheduled indices. */
  idxinv = hx_array_index_scheduled(k, sz, dsched, nsched, sched);

  /* check that allocation was successful. */
  if (!idx || !idxinv) {
    /* raise an exception. */
    raise("failed to allocate %d+%d indices", nidx, nsched);
    return NULL;
  }

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

  /* free the scheduled index array. */
  free(idxinv);

  /* return the indices. */
  return idx;
}

/* hx_array_index_printfn(): core function used by hx_array_print() to write
 * the contents of a multidimensional index to standard error.
 * @k: the array size.
 * @arr: the array of indices.
 * @s: the array variable name.
 */
void hx_array_index_printfn (int k, int *arr, const char *s) {
  /* declare a required variable. */
  int ki;

  /* print the variable name. */
  fprintf(stderr, "%s[%d] = (", s, k);

  /* loop over the array elements. */
  for (ki = 0; ki < k; ki++) {
    /* print the index value. */
    fprintf(stderr, "%d", arr[ki]);

    /* print commas to separate values. */
    if (ki < k - 1)
      fprintf(stderr, ", ");
  }

  /* end the line. */
  fprintf(stderr, ")\n");
  fflush(stderr);
}

