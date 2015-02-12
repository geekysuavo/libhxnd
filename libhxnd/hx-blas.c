
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

/* hx_no_trans(): function to retrieve an element from a matrix.
 *  - transpose: false
 *  - conjugate: false
 * see hx_blas_trans() for more details on arguments.
 */
void hx_no_trans (hx_array *A, hx_scalar *aij, int i, int j) {
  /* copy the array element into the scalar. */
  hx_data_copy(A->x + i + j * A->sz[0], aij->x, A->n);
}

/* hx_trans(): function to retrieve an element from a matrix.
 *  - transpose: true
 *  - conjugate: false
 * see hx_blas_trans() for more details on arguments.
 */
void hx_trans (hx_array *A, hx_scalar *aij, int i, int j) {
  /* copy the array element into the scalar. */
  hx_data_copy(A->x + j + i * A->sz[0], aij->x, A->n);
}

/* hx_conj_trans(): function to retrieve an element from a matrix.
 *  - transpose: true
 *  - conjugate: true
 * see hx_blas_trans() for more details on arguments.
 */
void hx_conj_trans (hx_array *A, hx_scalar *aij, int i, int j) {
  /* copy the array element into the scalar. */
  hx_data_conj(A->x + j + i * A->sz[0], aij->x, A->n);
}

