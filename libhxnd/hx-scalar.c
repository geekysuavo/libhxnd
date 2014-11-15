
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

/* hx_scalar_alloc(): allocate a hypercomplex scalar structure for a given
 * dimensionality. a pointer to the scalar structure is required by this
 * function.
 */
int hx_scalar_alloc (hx_scalar *x, int d) {
  /* check if the specified dimensionality is supported. */
  if (d < 0)
    throw("invalid algebraic dimensionality %d", d);

  /* store the dimensionality (d) and number of coefficients (n). */
  x->d = d;
  x->n = 1 << d;

  /* ensure that the d-dimensional shared multiplication table has been
   * initialized, and return failure if not.
   */
  if (!(x->tbl = hx_algebras_get(d)))
    throw("failed to retrieve %d-algebra", d);

  /* allocate the array of coefficients. return failure if allocation fails. */
  x->x = (real*) calloc(x->n, sizeof(real));
  if (x->x == NULL)
    throw("failed to allocate coefficient array");

  /* return success. */
  return 1;
}

/* hx_scalar_free(): de-allocate a previously allocated hypercomplex scalar
 * structure. a pointer to the scalar structure is required by this function.
 */
void hx_scalar_free (hx_scalar *x) {
  /* do not attempt to free a null pointer. */
  if (x == NULL)
    return;

  /* check if the coefficient array is allocated. */
  if (x->x != NULL) {
    /* yes. free its associated memory. */
    free(x->x);
    x->x = NULL;
  }

  /* de-initialize the dimensionality and coefficient count. */
  x->d = 0;
  x->n = 0;
}

/* hx_scalar_resize(): change the dimensionality of a hypercomplex scalar.
 * @x: a pointer to the scalar to resize.
 * @d: the new dimensionality.
 */
int hx_scalar_resize (hx_scalar *x, int d) {
  /* define a required variable. */
  int n;

  /* check if the specified dimensionality is supported. */
  if (d < 0)
    throw("invalid algebraic dimensionality %d", d);

  /* check that the new size is different. */
  if (d == x->d)
    return 1;

  /* compute the new number of coefficients. */
  n = 1 << d;

  /* reallocate the coefficients array in-place. */
  x->x = (real*) realloc(x->x, n * sizeof(real));

  /* check that the array was successfully reallocated. */
  if (!x->x)
    throw("failed to reallocate %d reals", n);

  /* initialize the new coefficients, if a grow was performed. */
  if (n > x->n)
    memset(x->x + x->n, 0, sizeof(real) * (n - x->n));

  /* store the new dimensionality constants. */
  x->d = d;
  x->n = n;

  /* ensure that the d-dimensional shared multiplication table has been
   * initialized, and return failure if not.
   */
  if (!(x->tbl = hx_algebras_get(d)))
    throw("failed to retrieve %d-algebra", d);

  /* return success. */
  return 1;
}

