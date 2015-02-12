
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

/* hx_blas_gemm(): compute the matrix-matrix product of two arrays.
 * @tA: whether to transpose the matrix @A.
 * @tB: whether to transpose the matrix @B.
 * @alpha: scale factor for the product.
 * @A: first matrix-shaped argument.
 * @B: second matrix-shaped argument.
 * @beta: second term scale factor.
 * @C: matrix shaped output array argument.
 */
int hx_blas_gemm (int tA, int tB, real alpha, hx_array *A, hx_array *B,
                  real beta, hx_array *C) {
  /* declare a few required variables:
   * @i, @j, @k: loop counters.
   * @n: number of array coefficients per scalar.
   * @M, @N: rows and columns of the output array.
   * @K: inner size (row/column) of the AB product.
   * @sum: temporary hypercomplex sum.
   */
  int i, j, k, n, M, N, K;
  int idxa, idxb, idxc;
  hx_scalar sum;

  /* ensure the arrays are all of correct shape. */
  hx_array_assert_matrix(A);
  hx_array_assert_matrix(B);
  hx_array_assert_matrix(C);

  /* ensure the arrays are all of correct dimensionality. */
  if (hx_array_dims_cmp(A, B) || hx_array_dims_cmp(A, C))
    throw("algebraic dimensionality mismatch");

  /* compute the output array sizes. */
  M = hx_matrix_rows(C);
  N = hx_matrix_cols(C);

  /* check if @C must be modified by beta. */
  if (beta == 0.0)
    hx_array_zero(C);
  else if (beta != 1.0)
    hx_blas_scal(beta, C);

  /* check if multiplication is even required. */
  if (alpha == 0.0)
    return 1;

  /* compute the number of coefficients per scalar. */
  n = A->n;

  /* check the transpose mode. */
  if (tA && tB) {
    /* get the inner loop index. */
    K = hx_matrix_cols(A);

    /* ensure the sizes match. */
    if (hx_matrix_rows(C) != hx_matrix_cols(A) ||
        hx_matrix_cols(C) != hx_matrix_rows(B) ||
        hx_matrix_rows(A) != hx_matrix_cols(B))
      throw("one or more operand size mismatches");
  }
  else if (tA) {
    /* get the inner loop index. */
    K = hx_matrix_rows(A);

    /* ensure the sizes match. */
    if (hx_matrix_rows(C) != hx_matrix_cols(A) ||
        hx_matrix_cols(C) != hx_matrix_cols(B) ||
        hx_matrix_rows(A) != hx_matrix_rows(B))
      throw("one or more operand size mismatches");
  }
  else if (tB) {
    /* get the inner loop index. */
    K = hx_matrix_cols(A);

    /* ensure the sizes match. */
    if (hx_matrix_rows(C) != hx_matrix_rows(A) ||
        hx_matrix_cols(C) != hx_matrix_rows(B) ||
        hx_matrix_cols(A) != hx_matrix_cols(B))
      throw("one or more operand size mismatches");
  }
  else {
    /* get the inner loop index. */
    K = hx_matrix_cols(A);

    /* ensure the sizes match. */
    if (hx_matrix_rows(C) != hx_matrix_rows(A) ||
        hx_matrix_cols(C) != hx_matrix_cols(B) ||
        hx_matrix_cols(A) != hx_matrix_rows(B))
      throw("one or more operand size mismatches");
  }

  /* allocate a scalar to hold intermediate sums. */
  if (!hx_scalar_alloc(&sum, A->d))
    throw("failed to allocate temporary %d-scalar", A->d);

  /* loop over the rows of the output matrix. */
  for (i = 0; i < M * n; i += n) {
    /* loop over the columns of the output matrix. */
    for (j = 0; j < N * n; j += n) {
      /* compute the output array linear index. */
      idxc = i + M * j;

      /* loop over the elements of the matrix product. */
      hx_scalar_zero(&sum);
      for (k = 0; k < K * n; k += n) {
        /* compute the first matrix linear index. */
        if (tA)
          idxa = i + A->sz[0] * k;
        else
          idxa = k + A->sz[0] * i;

        /* compute the second matrix linear index. */
        if (tB)
          idxb = k + B->sz[0] * j;
        else
          idxb = j + B->sz[0] * k;

        /* sum += A(i,k) * B(k,j) */
        hx_data_mul(A->x + idxa, B->x + idxb, sum.x, A->d, n, A->tbl);
      }

      /* scale the computed value and sum it into @C:
       *  C(i,j) += alpha * sum
       */
      hx_data_add(C->x + idxc, sum.x, C->x + idxc, alpha, C->d, C->n);
    }
  }

  /* free the temporary scalar. */
  hx_scalar_free(&sum);

  /* return success. */
  return 1;
}

