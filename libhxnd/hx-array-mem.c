
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

/* hx_array_init(): initialize the elements of a hypercomplex array.
 * @x: pointer to the array structure to initialize.
 */
void hx_array_init (hx_array *x) {
  /* initialize the structure members. */
  x->d = 0;
  x->n = 0;
  x->k = 0;
  x->len = 0;

  /* initialize the size array. */
  x->sz = NULL;

  /* initialize the coefficient array. */
  x->x = NULL;

  /* initialize the table pointer. */
  x->tbl = NULL;
}

/* hx_array_alloc(): allocate a hypercomplex array structure for a given
 * dimensionality. a pointer to the array structure is required by this
 * function.
 */
int hx_array_alloc (hx_array *x, int d, int k, hx_index sz) {
  /* declare a few required variables. */
  int i, len;

  /* check if the specified dimensionality is supported. */
  if (d < 0)
    throw("invalid algebraic dimensionality %d", d);

  /* check if the specified array dimensionality is supported. */
  if (k < 1)
    throw("invalid topological dimensionality %d", k);

  /* store the dimensionalities (d, k) and number of coefficients (n). */
  x->d = d;
  x->k = k;
  x->n = 1 << d;

  /* ensure that the d-dimensional shared multiplication table has been
   * initialized, and return failure if not.
   */
  if (!(x->tbl = hx_algebras_get(d)))
    throw("failed to retrieve %d-algebra", d);

  /* allocate the array of sizes. return failure if allocation fails. */
  x->sz = hx_index_alloc(x->k);
  if (x->sz == NULL)
    throw("failed to allocate size array");

  /* store the values in the size array locally. */
  for (i = 0, len = x->n; i < x->k; i++) {
    /* check that the size is valid. */
    if (sz[i] < 1)
      throw("dimension size %d (#%d) out of bounds [1,inf)", sz[i], i);

    /* store the size and multiply it into the total array length. */
    x->sz[i] = sz[i];
    len *= x->sz[i];
  }

  /* allocate the array of coefficients. fail if allocation fails. */
  x->x = (real*) calloc(len, sizeof(real));
  if (x->x == NULL)
    throw("failed to allocate coefficient array");

  /* store the array total coefficient count. */
  x->len = len;

  /* return success. */
  return 1;
}

/* hx_array_copy(): duplicate the contents of a hypercomplex array structure
 * into another structure.
 * @dst: destination array structure pointer.
 * @src: source array structure pointer.
 */
int hx_array_copy (hx_array *dst, hx_array *src) {
  /* allocate the destination array. */
  if (!hx_array_alloc(dst, src->d, src->k, src->sz))
    throw("failed to allocate destination array");

  /* copy the coefficients from the source array into the destination array.
   */
  memcpy(dst->x, src->x, src->len * sizeof(real));

  /* return success. */
  return 1;
}

/* hx_array_copy_real(): duplicate the real contents of a hypercomplex array
 * structure into another structure.
 * @dst: destination array structure pointer.
 * @src: source array structure pointer.
 */
int hx_array_copy_real (hx_array *dst, hx_array *src) {
  /* declare a few required variables:
   * @i: scalar array index.
   * @n: scalar array count.
   */
  int i, n;

  /* check if the source array is real. */
  if (hx_array_is_real(src))
    return hx_array_copy(dst, src);

  /* allocate the destination array. */
  if (!hx_array_alloc(dst, 0, src->k, src->sz))
    throw("failed to allocate destination array");

  /* compute the number of scalars in the source array. */
  n = dst->len;

  /* copy the real coefficients from the source array. */
  for (i = 0; i < n; i++)
    memcpy(dst->x + i, src->x + i * src->n, sizeof(real));

  /* return success. */
  return 1;
}

/* hx_array_free(): de-allocate a previously allocated hypercomplex array
 * structure. a pointer to the array structure is required by this function.
 */
void hx_array_free (hx_array *x) {
  /* do not attempt to free a null pointer. */
  if (x == NULL)
    return;

  /* free the coefficient array. */
  if (x->x)
    free(x->x);

  /* free the size array. */
  hx_index_free(x->sz);

  /* re-initialize the array structure contents. */
  hx_array_init(x);
}

