
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

/* ensure once-only inclusion. */
#ifndef __HXND_HX_BLAS_H__
#define __HXND_HX_BLAS_H__

/* hx_blas_trans(): function prototype for matrix transposition operations.
 * @A: pointer to the array for the input matrix.
 * @aij: pointer to the array for the output scalar.
 * @i: the row index of the final (transposed or not) matrix.
 * @j: the column index of the final (transposed or not) matrix.
 */
typedef void (*hx_blas_trans) (hx_array *A, hx_scalar *aij, int i, int j);

/* function declarations (hx-blas.c): */

void hx_no_trans (hx_array *A, hx_scalar *aij, int i, int j);

void hx_trans (hx_array *A, hx_scalar *aij, int i, int j);

void hx_conj_trans (hx_array *A, hx_scalar *aij, int i, int j);

/* function declarations (hx-blas-l1.c): */

int hx_blas_rdot (hx_array *x, hx_array *y, real *delta);

int hx_blas_cdotu (hx_array *x, hx_array *y, hx_scalar *delta);

int hx_blas_cdotc (hx_array *x, hx_array *y, hx_scalar *delta);

real hx_blas_nrm2 (hx_array *x);

real hx_blas_asum (hx_array *x);

int hx_blas_iamax (hx_array *x);

int hx_blas_swap (hx_array *x, hx_array *y);

int hx_blas_copy (hx_array *x, hx_array *y);

int hx_blas_scal (real alpha, hx_array *x);

int hx_blas_axpy (real alpha, hx_array *x, hx_array *y);

/* function declarations (hx-blas-l2.c): */

int hx_blas_gemv (hx_blas_trans tA, real alpha,
                  hx_array *A, hx_array *x,
                  real beta, hx_array *y);

int hx_blas_rger (real alpha, hx_array *x, hx_array *y, hx_array *A);

int hx_blas_cgeru (real alpha, hx_array *x, hx_array *y, hx_array *A);

int hx_blas_cgerc (real alpha, hx_array *x, hx_array *y, hx_array *A);

/* function declarations (hx-blas-l3.c): */

int hx_blas_gemm (hx_blas_trans tA, hx_blas_trans tB, real alpha,
                  hx_array *A, hx_array *B,
                  real beta, hx_array *C);

#endif /* __HXND_HX_BLAS_H__ */

