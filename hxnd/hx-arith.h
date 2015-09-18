
/* hxnd: A framework for n-dimensional hypercomplex calculations for NMR.
 * Copyright (C) 2014-2015  Bradley Worley  <geekysuavo@gmail.com>.
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
#ifndef __HXND_HX_ARITH_H__
#define __HXND_HX_ARITH_H__

/* function declarations, raw coefficient data: */

int hx_data_add (real *xa, real *xb, real *xc, real s, int d, int n);

int hx_data_mul (real *xa, real *xb, real *xc, int d, int n, hx_algebra tbl);

int hx_data_shuf (real *xa, real *xb, real *xc, real *xd,
                  real *xph, real *xtmp, int d, int n,
                  hx_algebra tbl);

int hx_data_copy (real *x, real *xcpy, int n);

int hx_data_conj (real *x, real *xh, int n);

int hx_data_zero (real *x, int n);

int hx_data_fill (real *x, int n, real val);

int hx_data_norm (real *x, int n);

real hx_data_real_norm (real *x, int n);

real hx_data_real_sumsq (real *x, int n);

int hx_data_negate_basis (real *x, int d, int n, int dneg);

int hx_data_reorder_bases (real *x, int d, int n, int *order);

/* function declarations, scalars: */

int hx_scalar_add (hx_scalar *a, hx_scalar *b, real s, hx_scalar *c);

int hx_scalar_mul (hx_scalar *a, hx_scalar *b, hx_scalar *c);

int hx_scalar_scale (hx_scalar *a, real s, hx_scalar *b);

int hx_scalar_zero (hx_scalar *a);

int hx_scalar_fill (hx_scalar *a, real val);

int hx_scalar_norm (hx_scalar *a);

int hx_scalar_negate_basis (hx_scalar *x, int dneg);

int hx_scalar_reorder_bases (hx_scalar *x, int *order);

/* function declarations, arrays: */

int hx_array_add_scalar (hx_array *a, hx_scalar *b, real s, hx_array *c);

int hx_array_add_array (hx_array *a, hx_array *b, real s, hx_array *c);

int hx_array_mul_scalar (hx_array *a, hx_scalar *b, hx_array *c);

int hx_array_mul_array (hx_array *a, hx_array *b, hx_array *c);

int hx_array_mul_vector (hx_array *a, hx_array *b, int kmul, hx_array *c);

int hx_array_scale (hx_array *a, real s, hx_array *b);

int hx_array_zero (hx_array *a);

int hx_array_fill (hx_array *a, real val);

int hx_array_norm (hx_array *a);

int hx_array_alternate_sign (hx_array *x, int k);

int hx_array_negate_basis (hx_array *x, int dneg);

int hx_array_reorder_bases (hx_array *x, int *order);

#endif /* __HXND_HX_ARITH_H__ */

