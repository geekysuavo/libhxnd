
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

/* * * * * * * * * * * LEVEL 1 * * * * * * * * * * */

/* hx_blas_dot(): compute the inner product of two arrays.
 * @x: first array operand.
 * @y: second array operand.
 */
real hx_blas_dot (hx_array *x, hx_array *y) {
  /* ensure the arrays are of equal length. */
  if (x->len != y->len)
    return 0.0;

  /* compute and return the result. */
#ifdef HX_DOUBLE_PRECISION
  return cblas_ddot(x->len, x->x, x->n, y->x, y->n);
#else
  return cblas_sdot(x->len, x->x, x->n, y->x, y->n);
#endif
}

/* hx_blas_nrm2(): compute the euclidean norm of an array.
 * @x: input array operand.
 */
real hx_blas_nrm2 (hx_array *x) {
  /* compute and return the result. */
#ifdef HX_DOUBLE_PRECISION
  return cblas_dnrm2(x->len, x->x, x->n);
#else
  return cblas_snrm2(x->len, x->x, x->n);
#endif
}

/* hx_blas_asum(): compute the sum of absolute values of an array.
 * @x: input array operand.
 */
real hx_blas_asum (hx_array *x) {
  /* compute and return the result. */
#ifdef HX_DOUBLE_PRECISION
  return cblas_dasum(x->len, x->x, x->n);
#else
  return cblas_sasum(x->len, x->x, x->n);
#endif
}

/* hx_blas_iamax(): locate the element of an array with the largest absolute
 * value.
 * @x: input array operand.
 */
int hx_blas_iamax (hx_array *x) {
  /* compute and return the result. */
#ifdef HX_DOUBLE_PRECISION
  return cblas_idamax(x->len, x->x, x->n);
#else
  return cblas_isamax(x->len, x->x, x->n);
#endif
}

/* hx_blas_swap(): swap the elements of two arrays.
 * @x: first array operand.
 * @y: second array operand.
 */
void hx_blas_swap (hx_array *x, hx_array *y) {
  /* ensure the arrays are of equal length. */
  if (x->len != y->len)
    return;

  /* execute the swap function. */
#ifdef HX_DOUBLE_PRECISION
  cblas_dswap(x->len, x->x, x->n, y->x, y->n);
#else
  cblas_sswap(x->len, x->x, x->n, y->x, y->n);
#endif
}

/* hx_blas_copy(): copy the elements of one array into another.
 * @x: first array operand.
 * @y: second array operand.
 */
void hx_blas_copy (hx_array *x, hx_array *y) {
  /* ensure the arrays are of equal length. */
  if (x->len != y->len)
    return;

  /* execute the copy function. */
#ifdef HX_DOUBLE_PRECISION
  cblas_dcopy(x->len, x->x, x->n, y->x, y->n);
#else
  cblas_scopy(x->len, x->x, x->n, y->x, y->n);
#endif
}

/* hx_blas_scal(): scale the elements of an array by a scalar factor.
 * @alpha: scale factor.
 * @x: array operand.
 */
void hx_blas_scal (real alpha, hx_array *x) {
  /* execute the scal function. */
#ifdef HX_DOUBLE_PRECISION
  cblas_dscal(x->len, alpha, x->x, x->n);
#else
  cblas_sscal(x->len, alpha, x->x, x->n);
#endif
}

/* hx_blas_axpy(): compute the sum or difference of two arrays.
 * @alpha: scale factor of first operand.
 * @x: first array operand.
 * @y: second array operand.
 */
void hx_blas_axpy (real alpha, hx_array *x, hx_array *y) {
  /* ensure the arrays are of equal length. */
  if (x->len != y->len)
    return;

  /* execute the axpy function. */
#ifdef HX_DOUBLE_PRECISION
  cblas_daxpy(x->len, alpha, x->x, x->n, y->x, y->n);
#else
  cblas_saxpy(x->len, alpha, x->x, x->n, y->x, y->n);
#endif
}

/* * * * * * * * * * * LEVEL 2 * * * * * * * * * * */

/* hx_blas_gemv(): compute the matrix-vector product of two arrays.
 * @tA: whether to transpose the matrix @A.
 * @A: matrix-shaped array operand.
 * @x: vector-shaped array operand.
 * @beta: second term scale factor.
 * @y: vector-shaped output array operand.
 */
void hx_blas_gemv (int tA, real alpha, hx_array *A, hx_array *x,
                   real beta, hx_array *y) {
  /* ensure the arrays are all of correct configuration. */
  if (A->k != 2 || x->k != 1 || y->k != 1)
    return;

  /* check the operand sizes. */
  if (tA) {
    /* with transposition. */
    if (A->sz[1] != y->sz[0] ||
        A->sz[0] != x->sz[0])
      return;
  }
  else {
    /* without transposition. */
    if (A->sz[0] != y->sz[0] ||
        A->sz[1] != x->sz[0])
      return;
  }

  /* execute the gemv function. */
#ifdef HX_DOUBLE_PRECISION
  cblas_dgemv(CblasColMajor, tA ? CblasTrans : CblasNoTrans,
              A->sz[0], A->sz[1], alpha, A->x, A->sz[0],
              x->x, x->n, beta, y->x, y->n);
#else
  cblas_sgemv(CblasColMajor, tA ? CblasTrans : CblasNoTrans,
              A->sz[0], A->sz[1], alpha, A->x, A->sz[0],
              x->x, x->n, beta, y->x, y->n);
#endif
}

/* hx_blas_ger(): compute the rank-1 update of an array by two arrays.
 * @alpha: scale factor for the update.
 * @x: first vector-shaped array operand.
 * @y: second vector-shaped array operand.
 * @A: matrix-shaped array operand.
 */
void hx_blas_ger (real alpha, hx_array *x, hx_array *y, hx_array *A) {
  /* ensure the arrays are all of correct configuration. */
  if (x->k != 1 || y->k != 1 || A->k != 2 ||
      A->sz[0] != x->sz[0] ||
      A->sz[1] != y->sz[0])
    return;

  /* execute the ger function. */
#ifdef HX_DOUBLE_PRECISION
  cblas_dger(CblasColMajor, A->sz[0], A->sz[1],
             alpha, x->x, x->n, y->x, y->n,
             A->x, A->sz[0]);
#else
  cblas_sger(CblasColMajor, A->sz[0], A->sz[1],
             alpha, x->x, x->n, y->x, y->n,
             A->x, A->sz[0]);
#endif
}

/* * * * * * * * * * * LEVEL 3 * * * * * * * * * * */

/* hx_blas_gemm(): compute the matrix-matrix product of two arrays.
 * @tA: whether to transpose the matrix @A.
 * @tB: whether to transpose the matrix @B.
 * @alpha: scale factor for the product.
 * @A: first matrix-shaped argument.
 * @B: second matrix-shaped argument.
 * @beta: second term scale factor.
 * @C: matrix shaped output array argument.
 */
void hx_blas_gemm (int tA, int tB, real alpha, hx_array *A, hx_array *B,
                   real beta, hx_array *C) {
  /* declare a few required variables. */
  int m, n, k;

  /* ensure the arrays are all of correct configuration. */
  if (A->k != 2 || B->k != 2 || C->k != 2)
    return;

  /* compute the array sizes. */
  m = B->sz[0];
  n = B->sz[1];
  k = (tA ? A->sz[0] : A->sz[1]);

  /* check the operand sizes. */
  if (tA && tB) {
    /* with transposition of both @A and @B. */
    if (A->sz[1] != C->sz[0] ||
        A->sz[0] != B->sz[1] ||
        B->sz[0] != C->sz[1])
      return;
  }
  else if (tA) {
    /* with transposition of @A. */
    if (A->sz[1] != C->sz[0] ||
        A->sz[0] != B->sz[0] ||
        B->sz[1] != C->sz[1])
      return;
  }
  else if (tB) {
    /* with transposition of @B. */
    if (A->sz[0] != C->sz[0] ||
        A->sz[1] != B->sz[1] ||
        B->sz[0] != C->sz[1])
      return;
  }
  else {
    /* without transposition. */
    if (A->sz[0] != C->sz[0] ||
        A->sz[1] != B->sz[0] ||
        B->sz[1] != C->sz[1])
      return;
  }

  /* execute the gemm function. */
#ifdef HX_DOUBLE_PRECISION
  cblas_dgemm(CblasColMajor,
              tA ? CblasTrans : CblasNoTrans,
              tB ? CblasTrans : CblasNoTrans,
              m, n, k, alpha, A->x, A->sz[0],
              B->x, B->sz[0], beta,
              C->x, C->sz[0]);
#else
  cblas_sgemm(CblasColMajor,
              tA ? CblasTrans : CblasNoTrans,
              tB ? CblasTrans : CblasNoTrans,
              m, n, k, alpha, A->x, A->sz[0],
              B->x, B->sz[0], beta,
              C->x, C->sz[0]);
#endif
}

