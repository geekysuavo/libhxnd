
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
#ifndef __HXND_HX_ARITH_H__
#define __HXND_HX_ARITH_H__

/* function declarations, raw coefficient data: */

int hx_data_add (real *xa, real *xb, real *xc, real s, int d, int n);

int hx_data_mul (real *xa, real *xb, real *xc, int d, int n, hx_algebra tbl);

int hx_data_norm (real *x, int d, int n);

/* function declarations, scalars: */

int hx_scalar_add (hx_scalar *a, hx_scalar *b, real s, hx_scalar *c);

int hx_scalar_mul (hx_scalar *a, hx_scalar *b, hx_scalar *c);

int hx_scalar_norm (hx_scalar *a);

/* function declarations, arrays: */

int hx_array_add_scalar (hx_array *a, hx_scalar *b, real s, hx_array *c);

int hx_array_add_array (hx_array *a, hx_array *b, real s, hx_array *c);

int hx_array_mul_scalar (hx_array *a, hx_scalar *b, hx_array *c);

int hx_array_mul_array (hx_array *a, hx_array *b, hx_array *c);

int hx_array_norm (hx_array *a);

#endif /* __HXND_HX_ARITH_H__ */

