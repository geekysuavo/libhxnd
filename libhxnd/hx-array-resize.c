
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

/* hx_array_is_real(): return whether a specified hypercomplex array contains
 * only real data.
 * @x: a pointer to the array to query.
 */
int hx_array_is_real (hx_array *x) {
  /* return whether the array is real or not. */
  return (x && x->d == 0);
}

/* hx_array_resize_d(): change the algebraic dimensionality of a
 * hypercomplex array. this function performs the resize operation
 * in-place.
 * @x: a pointer to the array to resize.
 * @d: the new algebraic dimensionality.
 */
int hx_array_resize_d (hx_array *x, int d) {
  /* declare a few required variables:
   * @is: scalar array index.
   * @ns: scalar array count.
   * @n: new number of coefficients.
   * @ncpy: number of bytes to copy per scalar.
   */
  int is, ns, n, nmin, ncpy;
  hx_scalar xs;

  /* check that the specified dimensionality is supported. */
  if (d < 0)
    throw("algebraic dimensionality %d is invalid", d);

  /* compute the new number of coefficients. */
  n = 1 << d;

  /* compute the new number of bytes per scalar. */
  nmin = (n < x->n ? n : x->n);
  ncpy = nmin * sizeof(real);

  /* compute the number of scalars in the array. */
  ns = x->len / x->n;

  /* allocate a temporary scalar. */
  if (!hx_scalar_alloc(&xs, d))
    throw("failed to allocate temporary %d-scalar", d);

  /* determine whether a shrink or a grow is required. */
  if (d < x->d) {
    /* shrink: loop forward through the scalar values of the array. */
    for (is = 0; is < ns; is++) {
      /* copy the scalar value from the old location to the new one. */
      memcpy(xs.x, x->x + is * x->n, ncpy);
      memcpy(x->x + is * n, xs.x, ncpy);
    }

    /* resize the coefficient array in place. */
    x->x = (real*) realloc(x->x, ns * n * sizeof(real));

    /* check that reallocation succeeded. */
    if (x->x == NULL)
      throw("failed to reallocate coefficient array");
  }
  else if (d > x->d) {
    /* resize the coefficient array in place. */
    x->x = (real*) realloc(x->x, ns * n * sizeof(real));

    /* check that reallocation succeeded. */
    if (x->x == NULL)
      throw("failed to reallocate coefficient array");

    /* grow: loop backward through the scalar values of the array. */
    for (is = ns - 1; is >= 0; is--) {
      /* copy the scalar value from the old location to the new one. */
      memcpy(xs.x, x->x + is * x->n, ncpy);
      memset(x->x + is * n, 0, n * sizeof(real));
      memcpy(x->x + is * n, xs.x, ncpy);
    }
  }

  /* store the new array dimensionality and length. */
  x->d = d;
  x->n = n;
  x->len = ns * n;

  /* free the temporary scalar. */
  hx_scalar_free(&xs);

  /* return success. */
  return 1;
}

/* hx_array_resize(): change the configuration of a hypercomplex array.
 * @x: a pointer to the array to resize.
 * @d: the new algebraic dimensionality.
 * @k: the new topological dimensionality.
 * @sz: the new size array.
 */
int hx_array_resize (hx_array *x, int d, int k, hx_index sz) {
  /* define a few required variables:
   * @arr: array of indices for iteration.
   * @idx: new linear iteration index.
   * @idxprev: old linear iteration index.
   * @i: general-purpose loop counter.
   * @n: new number of coefficients.
   * @len: new array length.
   * @ok: whether to copy a coefficient.
   * @nmin: smaller of new and old @n.
   * @kmax: larger of new and old @k.
   * @xnew: new coefficient array.
   */
  int pidx, pidxprev, i, n, len, ok;
  int nmin, kmax;
  hx_index idx;
  real *xnew;

  /* check if the array needs no resizing. */
  if (d == x->d && k == x->k && hx_index_cmp(k, sz, x->sz) == 0)
    return 1;

  /* check if only algebraic dimensionality is to be changed. */
  if (k == x->k && hx_index_cmp(k, sz, x->sz) == 0)
    return hx_array_resize_d(x, d);

  /* check if the specified dimensionalities are supported. */
  if (d < 0 || k < 1)
    throw("dimensionalities (%d, %d) are invalid", d, k);

  /* compute the new number of coefficients. */
  n = 1 << d;

  /* compute the sizes that are smaller, before or after. */
  nmin = (n < x->n ? n : x->n);

  /* allocate an array to hold the iteration indices. */
  kmax = (k > x->k ? k : x->k);
  idx = hx_index_alloc(kmax);

  /* initialize the linear indices. */
  pidx = pidxprev = 0;

  /* check that the index array was successfully allocated. */
  if (!idx)
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
  for (pidx = 0; pidx < len / n;) {
    /* loop over the dimensions to check if we're in bounds. */
    for (i = 0, ok = 1; i < k; i++) {
      /* check if we've exceeded one of the following:
       *  - the size of the currently indexed (new) dimension.
       *  - the dimensionality of the original array.
       *  - the size of the currently indexed (original) dimension.
       */
      if (idx[i] >= sz[i] ||
          (idx[i] && (i >= x->k || idx[i] >= x->sz[i]))) {
        /* yes. break the loop. */
        ok = 0;
        break;
      }
    }

    /* check if we're on a copiable index. */
    if (ok) {
      /* pack the indices into the old array linear index. */
      hx_index_pack(x->k, x->sz, idx, &pidxprev);

      /* loop over the coefficients. */
      for (i = 0; i < nmin; i++)
        xnew[i + n * pidx] = x->x[i + x->n * pidxprev];
    }

    /* increment the multidimensional indices. */
    hx_index_incr(k, sz, idx);
    pidx++;
  }

  /* free the allocated index array. */
  hx_index_free(idx);

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

