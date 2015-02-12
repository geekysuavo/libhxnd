
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

/* hx_blas_dot(): compute the inner product of two arrays.
 * @x: first array operand.
 * @y: second array operand.
 * @delta: output dot product scalar.
 */
int hx_blas_dot (hx_array *x, hx_array *y, hx_scalar *delta) {
  /* declare a required variable:
   * @i: array traversal loop index.
   */
  int i;

  /* ensure the arrays are of equal length. */
  if (x->len != y->len)
    throw("array length mismatch (%d != %d)", x->len, y->len);

  /* ensure the arrays are of equal dimensionality. */
  if (x->d != y->d || x->d != delta->d)
    throw("array dimensionality mismatch (%d != %d, %d)",
          x->d, y->d, delta->d);

  /* use a faster real-only function whenever possible. */
  if (hx_array_is_real(x)) {
    /* compute the real dot product. */
    for (i = 0, delta->x[0] = 0.0; i < x->len; i++)
      delta->x[0] += x->x[i] * y->x[i];
  }
  else {
    /* compute the hypercomplex dot product. */
    hx_scalar_zero(delta);
    for (i = 0; i < x->len; i += x->n)
      hx_data_mul(x->x + i, y->x + i, delta->x, x->d, x->n, x->tbl);
  }

  /* return success. */
  return 1;
}

/* hx_blas_sumsq(): compute the sum of squares of an array.
 * @x: input array operand.
 */
real hx_blas_sumsq (hx_array *x) {
  /* declare a few required variables:
   * @ssq: sum of squares of all array coefficients.
   * @i: array traversal loop index.
   */
  real ssq;
  int i;

  /* compute the array sum of squares. */
  for (i = 0, ssq = 0.0; i < x->len; i++)
    ssq += x->x[i] * x->x[i];

  /* return the computed result. */
  return ssq;
}

/* hx_blas_nrm2(): compute the euclidean norm of an array.
 * @x: input array operand.
 */
real hx_blas_nrm2 (hx_array *x) {
  /* return the square root of the array sum of squares. */
  return sqrt(hx_blas_sumsq(x));
}

/* hx_blas_asum(): compute the sum of absolute values of an array.
 * @x: input array operand.
 */
real hx_blas_asum (hx_array *x) {
  /* declare a few required variables:
   * @sum: sum of all array coefficients.
   * @i: array traversal loop index.
   */
  real sum;
  int i;

  /* compute the array sum of absolute values. */
  for (i = 0, sum = 0.0; i < x->len; i++)
    sum += fabs(x->x[i]);

  /* return the computed result. */
  return sum;
}

/* hx_blas_iamax(): locate the element of an array with the largest absolute
 * value.
 * @x: input array operand.
 */
int hx_blas_iamax (hx_array *x) {
  /* declare a few required variables:
   * @xi: currently indexed array value/result.
   * @xmax: maximum array value/result.
   * @i: array scalar element index.
   * @j: array coefficient index.
   * @imax: array index of @xmax.
   */
  int i, j, imax;
  real xi, xmax;

  /* loop over the scalar values of the array. */
  for (i = 0, imax = 0, xmax = 0.0; i < x->len; i += x->n) {
    /* compute the current array absolute value. */
    for (j = 0, xi = 0.0; j < x->n; j++)
      xi += fabs(x->x[i + j]);

    /* check if the current value exceeds the maximum value. */
    if (xi > xmax) {
      /* yes. store the current value and index. */
      xmax = xi;
      imax = i;
    }
  }

  /* return the identified index. */
  return imax / x->n;
}

/* hx_blas_swap(): swap the elements of two arrays.
 * @x: first array operand.
 * @y: second array operand.
 */
int hx_blas_swap (hx_array *x, hx_array *y) {
  /* declare a few required variables:
   * @swp: real scalar swap location.
   * @i: array traversal loop index.
   */
  real swp;
  int i;

  /* ensure the arrays are of equal length. */
  if (x->len != y->len)
    throw("array length mismatch (%d != %d)", x->len, y->len);

  /* loop over the array coefficients. */
  for (i = 0; i < x->len; i++) {
    /* swap the coefficients. */
    swp = x->x[i];
    x->x[i] = y->x[i];
    y->x[i] = swp;
  }

  /* return success. */
  return 1;
}

/* hx_blas_copy(): copy the elements of one array into another.
 * @x: source array operand.
 * @y: destination array operand.
 */
int hx_blas_copy (hx_array *x, hx_array *y) {
  /* ensure the arrays are of equal length. */
  if (x->len != y->len)
    throw("array length mismatch (%d != %d)", x->len, y->len);

  /* copy the array coefficients. */
  memcpy(y->x, x->x, x->len * sizeof(real));

  /* return success. */
  return 1;
}

/* hx_blas_scal(): scale the elements of an array by a scalar factor.
 * @alpha: scale factor.
 * @x: array operand.
 */
int hx_blas_scal (real alpha, hx_array *x) {
  /* declare a required variable:
   * @i: array traversal loop index.
   */
  int i;

  /* scale the coefficients of the array. */
  for (i = 0; i < x->len; i++)
    x->x[i] *= alpha;

  /* return success. */
  return 1;
}

/* hx_blas_axpy(): compute the sum or difference of two arrays.
 * @alpha: scale factor of first operand.
 * @x: first array operand.
 * @y: second array operand.
 */
int hx_blas_axpy (real alpha, hx_array *x, hx_array *y) {
  /* declare a required variable:
   * @i: array traversal loop index.
   */
  int i;

  /* ensure the arrays are of equal length. */
  if (x->len != y->len)
    throw("array length mismatch (%d != %d)", x->len, y->len);

  /* ensure the arrays are of equal dimensionality. */
  if (x->d != y->d)
    throw("array dimensionality mismatch (%d != %d)", x->d, y->d);

  /* compute the real vector sums. */
  for (i = 0; i < x->len; i++)
    y->x[i] += alpha * x->x[i];

  /* return success. */
  return 1;
}

