
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

/* hx_array_resize(): change the configuration of a hypercomplex array.
 * @x: a pointer to the array to resize.
 * @d: the new algebraic dimensionality.
 * @k: the new topological dimensionality.
 * @sz: the new size array.
 */
int hx_array_resize (hx_array *x, int d, int k, int *sz) {
  /* define a required variable. */
  int *arr, idx, idxprev, i, n, len, ok;
  int nmin, kmax;
  real *xnew;

  /* check if the specified dimensionalities are supported. */
  if (d < 0 || k < 1)
    throw("dimensionalities (%d, %d) are invalid", d, k);

  /* compute the new number of coefficients. */
  n = 1 << d;

  /* compute the sizes that are smaller, before or after. */
  nmin = (n < x->n ? n : x->n);

  /* allocate an array to hold the iteration indices. */
  kmax = (k > x->k ? k : x->k);
  arr = hx_array_index_alloc(kmax);
  idxprev = 0;

  /* check that the index array was successfully allocated. */
  if (!arr)
    throw("failed to allocate %d indices", kmax);

  /* compute the final size of the new coefficients array. */
  for (i = 0, len = n; i < k; i++) {
    /* check that the size is in bounds. */
    if (sz[i] < 1)
      throw("dimension size %d (#%d) out of bounds [1,inf)", sz[i], i);

    /* multiply the size into the total length. */
    len *= sz[i];
  }

  /* reallocate a brand new coefficients array. */
  xnew = (real*) calloc(len, sizeof(real));

  /* check that the data array was successfully allocated. */
  if (!xnew)
    throw("failed to allocate %d reals", len);

  /* loop over the set of indices. */
  for (idx = 0; idx < len / n;) {
    /* loop over the dimensions to check if we're in bounds. */
    for (i = 0, ok = 1; i < k; i++) {
      /* check if we've exceeded one of the following:
       *  - the size of the currently indexed (new) dimension.
       *  - the dimensionality of the original array.
       *  - the size of the currently indexed (original) dimension.
       */
      if (arr[i] >= sz[i] ||
          (arr[i] && (i >= x->k || arr[i] >= x->sz[i]))) {
        /* yes. break the loop. */
        ok = 0;
        break;
      }
    }

    /* check if we're on a copiable index. */
    if (ok) {
      /* pack the indices into the old array linear index. */
      hx_array_index_pack(x->k, x->sz, arr, &idxprev);

      /* loop over the coefficients. */
      for (i = 0; i < nmin; i++)
        xnew[i + n * idx] = x->x[i + x->n * idxprev];
    }

    /* increment the multidimensional indices. */
    hx_array_index_incr(k, sz, arr);
    idx++;
  }

  /* free the index arrays. */
  free(arr);

  /* store the new data array. */
  free(x->x);
  x->x = xnew;

  /* resize the size array, if the new size is different. */
  if (k != x->k)
    x->sz = (int*) realloc(x->sz, k * sizeof(int));

  /* store the new size array. */
  if (sz != x->sz)
    memcpy(x->sz, sz, k * sizeof(int));

  /* store the new dimensionality constants. */
  x->d = d;
  x->n = n;
  x->k = k;
  x->len = len;

  /* ensure that the d-dimensional shared multiplication table has been
   * initialized, and return failure if not.
   */
  if (!(x->tbl = hx_algebras_get(d)))
    throw("failed to retrieve %d-algebra", d);

  /* return success. */
  return 1;
}

