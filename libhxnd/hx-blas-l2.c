
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

/* hx_blas_gemv(): compute the matrix-vector product of two arrays.
 * @tA: whether to transpose the matrix @A.
 * @A: matrix-shaped array operand.
 * @x: vector-shaped array operand.
 * @beta: second term scale factor.
 * @y: vector-shaped output array operand.
 */
int hx_blas_gemv (int tA, real alpha, hx_array *A, hx_array *x,
                  real beta, hx_array *y) {
  /* declare a few required variables:
   * @i, @j: loop counters.
   * @idx: matrix linear index.
   * @sum: hypercomplex sum.
   */
  int i, j, idx;
  hx_scalar sum;

  /* ensure the arrays are all of correct shape. */
  hx_array_assert_matrix(A);
  hx_array_assert_vector(x);
  hx_array_assert_vector(y);

  /* ensure the arrays are all of correct dimensionality. */
  if (hx_array_dims_cmp(A, x) || hx_array_dims_cmp(A, y))
    throw("algebraic dimensionality mismatch");

  /* check if @y must be modified by beta. */
  if (beta == 0.0)
    hx_array_zero(y);
  else if (beta != 1.0)
    hx_blas_scal(beta, y);

  /* check if multiplication is even required. */
  if (alpha == 0.0)
    return 1;

  /* check the tranpose mode. */
  if (tA) {
    /* ensure the sizes match. */
    if (hx_vector_len(y) != hx_matrix_cols(A) ||
        hx_vector_len(x) != hx_matrix_rows(A))
      throw("one or more operand size mismatches");
  }
  else {
    /* ensure the sizes match. */
    if (hx_vector_len(y) != hx_matrix_rows(A) ||
        hx_vector_len(x) != hx_matrix_cols(A))
      throw("one or more operand size mismatches");
  }

  /* allocate a scalar to hold intermediate sums. */
  if (!hx_scalar_alloc(&sum, A->d))
    throw("failed to allocate temporary %d-scalar", A->d);

  /* loop over the elements of @y. */
  for (i = 0; i < y->len; i += y->n) {
    /* initialize the temporary sum. */
    hx_scalar_zero(&sum);

    /* loop over the elements of @x. */
    for (j = 0; j < x->len; j += x->n) {
      /* compute the matrix linear index. */
      if (tA)
        idx = j + A->sz[0] * i;
      else
        idx = i + A->sz[0] * j;

      /* sum += A(i,j) * x(i) */
      hx_data_mul(A->x + idx, x->x + j, sum.x, A->d, A->n, A->tbl);        
    }

    /* scale the computed value and sum it into @y:
     *  y(i) += alpha * sum
     */
    hx_data_add(y->x + i, sum.x, y->x + i, alpha, A->d, A->n);
  }

  /* free the temporary scalar. */
  hx_scalar_free(&sum);

  /* return success. */
  return 1;
}

/* hx_blas_ger(): compute the rank-1 update of an array by two arrays.
 * @alpha: scale factor for the update.
 * @x: first vector-shaped array operand.
 * @y: second vector-shaped array operand.
 * @A: matrix-shaped array operand.
 */
int hx_blas_ger (real alpha, hx_array *x, hx_array *y, hx_array *A) {
  /* declare a few required variables:
   * @i, @j: loop counters.
   * @idx: matrix linear index.
   * @hprod: hypercomplex temporary, unscaled product.
   */
  hx_scalar hprod;
  int i, j, idx;

  /* ensure the arrays are all of correct shape. */
  hx_array_assert_vector(x);
  hx_array_assert_vector(y);
  hx_array_assert_matrix(A);

  /* ensure the arrays are all of correct dimensionality. */
  if (hx_array_dims_cmp(x, A) || hx_array_dims_cmp(y, A))
    throw("algebraic dimensionality mismatch");

  /* ensure the array sizes match. */
  if (hx_vector_len(x) != hx_matrix_rows(A) ||
      hx_vector_len(y) != hx_matrix_cols(A))
    throw("one or more operand size mismatches");

  /* use a faster real-only matrix update whenever possible. */
  if (hx_array_is_real(A)) {
    /* update the elements of the real matrix. */
    for (j = idx = 0; j < y->len; j++)
      for (i = 0; i < x->len; i++, idx++)
        A->x[idx] += alpha * x->x[i] * y->x[j];
  }
  else {
    /* allocate a temporary hypercomplex scalar. */
    if (!hx_scalar_alloc(&hprod, A->d))
      throw("failed to allocate temporary %d-scalar", A->d);

    /* update the elements of the hypercomplex matrix. */
    for (j = idx = 0; j < y->len; j += y->n) {
      for (i = 0; i < x->len; i += x->n, idx += A->n) {
        /* compute the scaled product of the two vector entries,
         * then sum it with the appropriate matrix element.
         */
        hx_scalar_zero(&hprod);
        hx_data_mul(x->x + i, y->x + j, hprod.x, x->d, x->n, x->tbl);
        hx_data_add(A->x + idx, hprod.x, A->x + idx, alpha, A->d, A->n);
      }
    }

    /* free the temporary scalar. */
    hx_scalar_free(&hprod);
  }

  /* return success. */
  return 1;
}

