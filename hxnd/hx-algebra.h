
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

/* ensure once-only inclusion. */
#ifndef __HXND_HX_ALGEBRA_H__
#define __HXND_HX_ALGEBRA_H__

/* hx_algebra: multiplication table for hypercomplex nD scalars.
 *
 * this table is a row-major square array of signed ints that convey the
 * bases and signs of the products of all possible bases of a hypercomplex
 * number algebra.
 *
 * a d-dimensional hypercomplex table will contain n**2 elements, and will be
 * accessed in the following way:
 *
 * tbl[i * n + j] => sgn[i,j] * (k + 1) => sgn{ui * uj} * idx{uk}
 *
 * where k (idx{uk}) is the index of the output coefficient of the product,
 * and sgn[i,j] (sgn{ui * uj}) indicates whether the coefficient needs to
 * be sign-flipped (e.g. when i == j).
 */
typedef int *hx_algebra;

/* function declarations: */

void hx_algebras_init (void);

int hx_algebras_add (int d);

hx_algebra hx_algebras_get (int d);

#endif /* __HXND_HX_ALGEBRA_H__ */

