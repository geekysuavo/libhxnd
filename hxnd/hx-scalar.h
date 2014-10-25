
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
#ifndef __HXND_HX_SCALAR_H__
#define __HXND_HX_SCALAR_H__

/* hx_scalar: data type for hypercomplex nD scalars.
 *
 * a d-dimensional hypercomplex number is described by 2**d real coefficients,
 * where the order of the coefficients in 'x' is based on a binary indexing
 * system. for example:
 *
 * d = 0: x = {a},
 *  where x[0b] => a
 *
 * d = 1: x = {a + b u1},
 *  where x[0b] => a
 *        x[1b] => b
 *
 * d = 2: x = {a + b u1 + c u2 + d u1u2},
 *  where x[00b] => a
 *        x[01b] => b
 *        x[10b] => c
 *        x[11b] => d
 *
 * d = 3: x = {a + b u1 + c u2 + d u3 + e u1u2 + f u1u3 + g u2u3 + h u1u2u3},
 *  where x[000b] => a
 *        x[001b] => b
 *        x[010b] => c
 *        x[100b] => d
 *        x[011b] => e
 *        x[101b] => f
 *        x[110b] => g
 *        x[111b] => h
 *
 * in summary, a high-bit in a given position 'i' indicates that the 'i'-th
 * coefficient is multiplied by ui. it is important to notice that the real
 * coefficients WILL NOT be stored alphabetically in 'x'. instead, the first
 * 'd' coefficients are stored as 0, 1<<0, 1<<1, ... 1<<(d-1). indices of all
 * remaining basis products (u1u2, u1u3, etc) are then simple binary OR's of
 * the first 'd' bases (e.g. x{u1u3}: x[001b | 100b] => x[101b] => f).
 */
typedef struct {
  /* d: dimensionality.
   * n: number of coefficients (2**d).
   */
  int d, n;

  /* real coefficients. */
  real *x;

  /* multiplication table. this is a shared table, so it should not be freed
   * with the number.
   */
  hx_algebra tbl;
}
hx_scalar;

/* function declarations: */

int hx_scalar_alloc (hx_scalar *x, int d);

int hx_scalar_free (hx_scalar *x);

int hx_scalar_resize (hx_scalar *x, int d);

int hx_scalar_zero (hx_scalar *x);

int hx_scalar_phasor (hx_scalar *x, int d, real phi);

#endif /* __HXND_HX_SCALAR_H__ */

