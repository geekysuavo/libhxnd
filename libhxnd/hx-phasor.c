
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

/* hx_scalar_phasor(): stores the values of a phasor in a hypercomplex scalar.
 * @x: a pointer to the scalar to modify.
 * @d: the dimension index to act upon.
 * @phi: the phase angle, in [0, 2 pi].
 */
int hx_scalar_phasor (hx_scalar *x, int d, real phi) {
  /* declare a required variable. */
  int n;

  /* check that the dimension index is valid. */
  if (d < 0 || d >= x->d)
    throw("algebraic dimension %d out of bounds [0,%d)", d, x->d);

  /* compute the coefficient index for phasing. */
  n = 1 << d;

  /* zero the values of the coefficient array. */
  hx_scalar_zero(x);

  /* store the real and imaginary phase factors. */
  x->x[0] = cos(phi);
  x->x[n] = sin(phi);

  /* return success. */
  return 1;
}

/* hx_array_phasor(): stores the values of a phasor in a one-dimensional
 * hypercomplex array.
 * @x: a pointer to the array to modify.
 * @d: the dimension index to act upon.
 * @phi0: zero-order phase angle, in [0, 2 pi].
 * @phi1: first-order phase angle, unbounded.
 * @pivot: fractional pivot index, in [0, 1].
 */
int hx_array_phasor (hx_array *x, int d, real phi0, real phi1, real pivot) {
  /* declare a few required variables. */
  real fi, phi;
  int idx, n;

  /* check that the dimension index is valid. */
  if (d < 0 || d >= x->d)
    throw("algebraic dimension %d out of bounds [0,%d)", d, x->d);

  /* compute the coefficient index for phasing. */
  n = 1 << d;

  /* zero the values of the coefficient array. */
  hx_array_zero(x);

  /* loop over the array coefficients. */
  for (idx = 0; idx < x->len; idx += x->n) {
    /* compute the fractional array index. */
    fi = ((real) idx) / ((real) (x->len - 1));

    /* compute the current phase angle. */
    phi = phi0 + phi1 * (fi - pivot);

    /* store the current real and imaginary phase factors. */
    x->x[idx] = cos(phi);
    x->x[idx + n] = sin(phi);
  }

  /* return success. */
  return 1;
}

