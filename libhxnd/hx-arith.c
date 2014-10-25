
/* ndmath: A framework for n-dimensional hypercomplex calculations for NMR.
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

/* hx_data_add(): add the raw array elements of two hypercomplex values.
 * @xa: the raw array data of the first operand.
 * @xb: the raw array data of the second operand.
 * @xc: the raw array data of the result.
 * @s: the real factor to apply to b during addition.
 * @d: the algebraic dimensionality of the scalars.
 * @n: the number of array elements of the scalars.
 *
 * operation:
 *   c <= a + s * b
 *
 * operands:
 *   a: hypercomplex scalar.
 *   b: hypercomplex scalar.
 *   s: real scalar.
 *   c: hypercomplex scalar.
 */
int hx_data_add (real *xa, real *xb, real *xc, real s, int d, int n) {
  /* declare a required variable. */
  int i;

  /* check if a modified operation is to be performed. */
  if (xa == NULL) {
    /* modified: only scale b[i]. */
    for (i = 0; i < n; i++)
      xc[i] = s * xb[i];
  }
  else if (xb == NULL) {
    /* modified: only add s to a[i]. (i.e. b[i] == 1) */
    for (i = 0; i < n; i++)
      xc[i] = xa[i] + s;
  }
  else {
    /* full: scale the value of b[i] and add the result to a[i] */
    for (i = 0; i < n; i++)
      xc[i] = xa[i] + s * xb[i];
  }

  /* return success. */
  return 1;
}

/* hx_data_mul(): multiply the raw array elements of two hypercomplex values.
 * @xa: the raw array data of the first operand.
 * @xb: the raw array data of the second operand.
 * @xc: the raw array data of the result.
 * @d: the algebraic dimensionality of the scalars.
 * @n: the number of array elements of the scalars.
 * @tbl: the multiplication table to utilize.
 *
 * operation:
 *   c <= c + a * b
 *
 * operands:
 *   a: hypercomplex scalar.
 *   b: hypercomplex scalar.
 *   c: hypercomplex scalar.
 */
int hx_data_mul (real *xa, real *xb, real *xc, int d, int n, hx_algebra tbl) {
  /* declare a few required variables. */
  int i, j, k, tij;
  real sgn;

  /* loop over the first value's array elements. */
  for (i = 0; i < n; i++) {
    /* skip zero-valued left-operand coefficients. */
    if (xa[i] == 0.0)
      continue;

    /* loop over the second value's array elements. */
    for (j = 0; j < n; j++) {
      /* get the output index from the algebra table. */
      tij = tbl[i * n + j];
      k = (tij > 0 ? tij : -tij);

      /* get the output sign from the algebra table. */
      sgn = (real) (tij / k);
      k--;

      /* compute the final result. */
      xc[k] += sgn * xa[i] * xb[j];
    }
  }

  /* return success. */
  return 1;
}

/* hx_scalar_add(): add two hypercomplex scalar values.
 * @a: the structure pointer of the first operand.
 * @b: the structure pointer of the second operand.
 * @s: the real factor to apply to b during addition.
 * @c: the structure pointer of the result.
 *
 * operation:
 *   c <= a + s * b
 *
 * operands:
 *   a: hypercomplex scalar.
 *   b: hypercomplex scalar.
 *   s: real scalar.
 *   c: hypercomplex scalar.
 */
int hx_scalar_add (hx_scalar *a, hx_scalar *b, real s, hx_scalar *c) {
  /* check if the algebraic dimensionalities match. */
  if (hx_scalar_dims_cmp(a, b) != 0 ||
      hx_scalar_dims_cmp(a, c) != 0)
    return 0;

  /* perform the raw data operation. */
  return hx_data_add(a->x, b->x, c->x, s, a->d, a->n);
}

/* hx_scalar_mul(): multiply two hypercomplex scalar values.
 * @a: the structure pointer of the first operand.
 * @b: the structure pointer of the second operand.
 * @c: the structure pointer of the result.
 *
 * operation:
 *   c <= c + a * b
 *
 * operands:
 *   a: hypercomplex scalar.
 *   b: hypercomplex scalar.
 *   c: hypercomplex scalar.
 */
int hx_scalar_mul (hx_scalar *a, hx_scalar *b, hx_scalar *c) {
  /* check if the algebraic dimensionalities match. */
  if (hx_scalar_dims_cmp(a, b) != 0 ||
      hx_scalar_dims_cmp(a, c) != 0)
    return 0;

  /* perform the raw data operation. */
  return hx_data_mul(a->x, b->x, c->x, a->d, a->n, a->tbl);
}

/* hx_array_add_scalar(): add a hypercomplex array and a scalar.
 * @a: the structure pointer to the array operand.
 * @b: the structure pointer to the scalar operand.
 * @s: the real factor to apply to b during addition.
 * @c: the structure pointer to the result.
 *
 * operation:
 *   c <= a + s * b
 *
 * operands:
 *   a: hypercomplex array.
 *   b: hypercomplex scalar.
 *   s: real scalar.
 *   c: hypercomplex array.
 */
int hx_array_add_scalar (hx_array *a, hx_scalar *b, real s, hx_array *c) {
  /* declare a required variable. */
  int i;

  /* check if the algebraic dimensionalities match. */
  if (a->d != b->d || hx_array_conf_cmp(a, c) != 0)
    return 0;

  /* loop over the array elements. */
  for (i = 0; i < a->len; i++) {
    /* perform the raw scalar data operation. */
    if (!hx_data_add(a->x + i, b->x, c->x + i, s, a->d, a->n))
      return 0;
  }

  /* return success. */
  return 1;
}

/* hx_array_add_array(): add two hypercomplex arrays.
 * @a: the structure pointer to the first operand.
 * @b: the structure pointer to the second operand.
 * @s: the real factor to apply to b during addition.
 * @c: the structure pointer to the result.
 *
 * operation:
 *   c <= a + s * b
 *
 * operands:
 *   a: hypercomplex array.
 *   b: hypercomplex array.
 *   s: real scalar.
 *   c: hypercomplex array.
 */
int hx_array_add_array (hx_array *a, hx_array *b, real s, hx_array *c) {
  /* declare a required variable. */
  int i;

  /* check if the array configurations match. */
  if (hx_array_conf_cmp(a, b) != 0 ||
      hx_array_conf_cmp(a, c) != 0)
    return 0;

  /* loop over the array elements. */
  for (i = 0; i < a->len; i += a->n) {
    /* perform the raw scalar data operation. */
    if (!hx_data_add(a->x + i, b->x + i, c->x +i, s, a->d, a->n))
      return 0;
  }

  /* return success. */
  return 1;
}

/* hx_array_mul_scalar(): multiply a hypercomplex array and a scalar.
 * @a: the structure pointer to the array operand.
 * @b: the structure pointer to the scalar operand.
 * @c: the structure pointer to the result.
 *
 * operation:
 *   c <= c + a * b
 *
 * operands:
 *   a: hypercomplex array.
 *   b: hypercomplex scalar.
 *   c: hypercomplex array.
 */
int hx_array_mul_scalar (hx_array *a, hx_scalar *b, hx_array *c) {
  /* declare a required variable. */
  int i;

  /* check if the algebraic dimensionalities match. */
  if (a->d != b->d || hx_array_conf_cmp(a, c) != 0)
    return 0;

  /* loop over the array elements. */
  for (i = 0; i < a->len; i++) {
    /* perform the raw scalar data operation. */
    if (!hx_data_mul(b->x, a->x + i, c->x + i, a->d, a->n, a->tbl))
      return 0;
  }

  /* return success. */
  return 1;
}

/* hx_array_mul_array(): multiply two hypercomplex arrays.
 * @a: the structure pointer to the first operand.
 * @b: the structure pointer to the second operand.
 * @c: the structure pointer to the result.
 *
 * operation:
 *   c <= c + a * b
 *
 * operands:
 *   a: hypercomplex array.
 *   b: hypercomplex array.
 *   c: hypercomplex array.
 */
int hx_array_mul_array (hx_array *a, hx_array *b, hx_array *c) {
  /* declare a required variable. */
  int i;

  /* check if the array configurations match. */
  if (hx_array_conf_cmp(a, b) != 0 ||
      hx_array_conf_cmp(a, c) != 0)
    return 0;

  /* loop over the array elements. */
  for (i = 0; i < a->len; i += a->n) {
    /* perform the raw scalar data operation. */
    if (!hx_data_mul(a->x + i, b->x + i, c->x +i, a->d, a->n, a->tbl))
      return 0;
  }

  /* return success. */
  return 1;
}

