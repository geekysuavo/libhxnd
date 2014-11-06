
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

/* algebras: centralized array of all defined hypercomplex algebras. once an
 * algebra of a certain dimensionality has been defined through a call to
 * hx_scalar_alloc() or hx_array_alloc(), it is simply referenced again
 * in subsequent calls that use the same dimensionality.
 */
hx_algebra *algebras;
int n_algebras;

/* hx_algebras_init(): initializes the shared array of hypercomplex
 * algebras. this must be called prior to handling any hypercomplex
 * numbers.
 */
void hx_algebras_init (void) {
  /* initialize the array to empty. */
  algebras = NULL;
  n_algebras = 0;
}

/* hx_algebras_add(): initializes a multiplication table having the specified
 * dimensionality (d) in the shared array of hypercomplex algebras. this must
 * be called prior to performing arithmetic on d-dimensional hypercomplex
 * numbers, but a manual call is not required: allocating a d-dimensional
 * scalar or array will automatically call this function.
 */
int hx_algebras_add (int d) {
  /* declare a few required variables. */
  int tij, i, j, k, n, nn, n_algebras_prev;

  /* check if the algebras array is large enough. */
  if (n_algebras < d + 1) {
    /* not large enough. we need to expand it. */
    n_algebras_prev = n_algebras;
    n_algebras = d + 1;

    /* reallocate the algebras array to the new size. */
    algebras = (hx_algebra*) realloc(algebras,
      n_algebras * sizeof(hx_algebra));

    /* check if the reallocation failed, and return failure if so. */
    if (algebras == NULL)
      throw("failed to resize algebras array");

    /* loop over the newly created algebras in the resized array. */
    for (i = n_algebras_prev; i < n_algebras; i++) {
      /* initialize the new algebra to null. */
      algebras[i] = NULL;
    }
  }

  /* check if the algebras array has added the desired dimensionality. */
  if (algebras[d] == NULL) {
    /* compute the size of the multiplication table array. */
    n = 1 << d;
    nn = n * n;

    /* allocate the multiplication table array. */
    algebras[d] = (hx_algebra) calloc(nn, sizeof(int));

    /* check if the allocation failed, and return failure if so. */
    if (algebras[d] == NULL)
      throw("failed to allocate %d-d algebra", d);

    /* loop over the first dimension (rows) of the multiplication table. */
    for (i = 0; i < n; i++) {
      /* loop over the second dimension (columns) of the table. */
      for (j = 0; j < n; j++) {
        /* compute the unsigned output index of the result. */
        tij = (i ^ j) + 1;

        /* check if the two complex basis elements will overlap to yield one
         * or more real negative factors (i.e. i*i=-1, j*j=-1, etc.).
         */
        for (k = 0; k < 8 * sizeof(int); k++) {
          /* check if the basis elements overlap in the currently indexed
           * dimension. if they overlap, negate the result.
           */
          if ((i & j) & (1 << k))
            tij *= -1;
        }

        /* store the signed output index of the result. */
        algebras[d][i * n + j] = tij;
      }
    }
  }

  /* return success. */
  return 1;
}

/* hx_algebras_get(): returns the requested algebras array, and builds it
 * if it hasn't been built already.
 */
hx_algebra hx_algebras_get (int d) {
  /* ensure that the d-dimensional shared multiplication table has been
   * initialized, and return failure if not.
   */
  if (!hx_algebras_add(d))
    return NULL;

  /* store a pointer to the d-dimensional multiplication table. */
  return algebras[d];
}

