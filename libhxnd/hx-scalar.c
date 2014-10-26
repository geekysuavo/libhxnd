
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
    return 0;

  /* store the dimensionality (d) and number of coefficients (n). */
  x->d = d;
  x->n = 1 << d;

  /* ensure that the d-dimensional shared multiplication table has been
   * initialized, and return failure if not.
   */
  if (!(x->tbl = hx_algebras_get(d)))
    return 0;

  /* allocate the array of coefficients. return failure if allocation fails. */
  x->x = (real*) calloc(x->n, sizeof(real));
  if (x->x == NULL)
    return 0;

  /* return success. */
  return 1;
}

/* hx_scalar_free(): de-allocate a previously allocated hypercomplex scalar
 * structure. a pointer to the scalar structure is required by this function.
 */
int hx_scalar_free (hx_scalar *x) {
  /* do not attempt to free a null pointer. */
  if (x == NULL)
    return 0;

  /* check if the coefficient array is allocated. */
  if (x->x != NULL) {
    /* yes. free its associated memory. */
    free(x->x);
    x->x = NULL;
  }

  /* de-initialize the dimensionality and coefficient count. */
  x->d = 0;
  x->n = 0;

  /* return success. */
  return 1;
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
    return 0;

  /* check that the new size is different. */
  if (d == x->d)
    return 1;

  /* compute the new number of coefficients. */
  n = 1 << d;

  /* reallocate the coefficients array in-place. */
  x->x = (real*) realloc(x->x, n * sizeof(real));

  /* check that the array was successfully reallocated. */
  if (!x->x)
    return 0;

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
    return 0;

  /* return success. */
  return 1;
}

/* hx_scalar_zero(): sets all coefficients in a scalar to zero.
 * @x: a pointer to the scalar to modify.
 */
int hx_scalar_zero (hx_scalar *x) {
  /* zero the scalar coefficient data as quickly as possible. */
  memset(x->x, 0, x->n * sizeof(real));

  /* return success. */
  return 1;
}

/* hx_scalar_phasor(): stores the values of a phasor in a hypercomplex scalar.
 * @x: a pointer to the scalar to modify.
 * @d: the dimension index to act upon.
 * @phi: the phase angle, in [0, 2 pi].
 */
int hx_scalar_phasor (hx_scalar *x, int d, real phi) {
  /* declare a required variable. */
  int n;

  /* check that the dimension index is valid. */
  if (d < 0 || d >= x->d)
    return 0;

  /* compute the coefficient index for phasing. */
  n = 1 << d;

  /* zero the values of the coefficient array. */
  memset(x->x, 0, x->n * sizeof(real));

  /* store the 'real' and 'imaginary' phase factors. */
  x->x[0] = cos(phi);
  x->x[n] = sin(phi);

  /* return success. */
  return 1;
}

