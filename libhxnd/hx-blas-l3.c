
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
 * @tA: transposition mode of the matrix @A.
 * @tB: transposition mode of the matrix @B.
 * @alpha: scale factor for the product.
 * @A: first matrix-shaped argument.
 * @B: second matrix-shaped argument.
 * @beta: second term scale factor.
 * @C: matrix shaped output array argument.
 */
int hx_blas_gemm (hx_blas_trans tA, hx_blas_trans tB, real alpha,
                  hx_array *A, hx_array *B,
                  real beta, hx_array *C) {
  /* declare a few required variables:
   * @i, @j, @k: loop counters.
   * @n: number of array coefficients per scalar.
   * @M, @N: rows and columns of the output array.
   * @K: inner size (row/column) of the AB product.
   * @sum: temporary hypercomplex sum.
   */
  int i, j, k, n, M, N, K, idx;
  hx_scalar a, b, sum;

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
  if (tA == &hx_no_trans && tB == &hx_no_trans) {
    /* get the inner loop index. */
    K = hx_matrix_cols(A);

    /* ensure the sizes match. */
    if (hx_matrix_rows(C) != hx_matrix_rows(A) ||
        hx_matrix_cols(C) != hx_matrix_cols(B) ||
        hx_matrix_cols(A) != hx_matrix_rows(B))
      throw("one or more operand size mismatches");
  }
  else if (tA == &hx_no_trans) {
    /* get the inner loop index. */
    K = hx_matrix_cols(A);

    /* ensure the sizes match. */
    if (hx_matrix_rows(C) != hx_matrix_rows(A) ||
        hx_matrix_cols(C) != hx_matrix_rows(B) ||
        hx_matrix_cols(A) != hx_matrix_cols(B))
      throw("one or more operand size mismatches");
  }
  else if (tB == &hx_no_trans) {
    /* get the inner loop index. */
    K = hx_matrix_rows(A);

    /* ensure the sizes match. */
    if (hx_matrix_rows(C) != hx_matrix_cols(A) ||
        hx_matrix_cols(C) != hx_matrix_cols(B) ||
        hx_matrix_rows(A) != hx_matrix_rows(B))
      throw("one or more operand size mismatches");
  }
  else {
    /* get the inner loop index. */
    K = hx_matrix_cols(A);

    /* ensure the sizes match. */
    if (hx_matrix_rows(C) != hx_matrix_cols(A) ||
        hx_matrix_cols(C) != hx_matrix_rows(B) ||
        hx_matrix_rows(A) != hx_matrix_cols(B))
      throw("one or more operand size mismatches");
  }

  /* allocate a scalar to hold intermediate sums. */
  if (!hx_scalar_alloc(&sum, A->d) ||
      !hx_scalar_alloc(&a, A->d) ||
      !hx_scalar_alloc(&b, A->d))
    throw("failed to allocate temporary %d-scalar", A->d);

  /* loop over the columns of the output matrix. */
  for (j = 0, idx = 0; j < N * n; j += n) {
    /* loop over the rows of the output matrix. */
    for (i = 0; i < M * n; i += n, idx += n) {
      /* loop over the elements of the matrix product. */
      hx_scalar_zero(&sum);
      for (k = 0; k < K * n; k += n) {
        /* extract the input matrix elements. */
        tA(A, &a, i, k);
        tB(B, &b, k, j);

        /* sum += A(i,k) * B(k,j) */
        hx_data_mul(a.x, b.x, sum.x, A->d, n, A->tbl);
      }

      /* scale the computed value and sum it into @C:
       *  C(i,j) += alpha * sum
       */
      hx_data_add(C->x + idx, sum.x, C->x + idx, alpha, C->d, C->n);
    }
  }

  /* free the temporary scalars. */
  hx_scalar_free(&sum);
  hx_scalar_free(&a);
  hx_scalar_free(&b);

  /* return success. */
  return 1;
}

