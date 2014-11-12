
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

/* hx_ispow2(): determines whether an integer is a power of two.
 * @value: integer to check.
 */
int hx_ispow2 (unsigned int value) {
  /* declare a few required variables:
   * @i: loop counter.
   * @n: high bit counter.
   */
  unsigned int i, n;

  /* just fail right away for the one wrong case. */
  if (value == 1)
    return 0;

  /* count the number of set bits in the value. */
  for (i = 0, n = 0; i < 8 * sizeof(value); i++)
    n += (value & (1 << i) ? 1 : 0);

  /* return true if the number of set bits is one. */
  return (n == 1);
}

/* hx_prevpow2(): finds the greatest power of two that is less than @value.
 * @value: the integer in question.
 */
unsigned int hx_prevpow2 (unsigned int value) {
  /* declare a few required variables. */
  unsigned int pow2;

  /* loop until we've found the value. */
  for (pow2 = 1; pow2 < value;)
    pow2 <<= 1;

  /* return the identified value. */
  pow2 >>= 1;
  return pow2;
}

/* hx_nextpow2(): finds the smallest power of two that is greater than @value.
 * @value: the integer in question.
 */
unsigned int hx_nextpow2 (unsigned int value) {
  /* declare a few required variables. */
  unsigned int pow2;

  /* loop until we've found the value. */
  for (pow2 = 1; pow2 <= value;)
    pow2 <<= 1;

  /* return the identified value. */
  return pow2;
}

/* hx_array_fft1d(): computes an in-place radix-2 fast fourier transform
 * of a hypercomplex one-dimensional array.
 *
 * args:
 *  see hx_array_vector_cb().
 *
 * varargs:
 *  @d: dimension to transform.
 *  @dir: direction of transformation.
 *  @w: preallocated hypercomplex scalar.
 *  @swp: preallocated hypercomplex scalar.
 */
int hx_array_fft1d (hx_array *x, hx_array *y,
                    int *arr, int idx,
                    va_list *vl) {
  /* declare a few required variables:
   * @i, @j, @k, @m: loop counters.
   * @n: array (scalar) element count.
   * @ncpy: number of bytes per scalar.
   * @step: loop stride.
   * @phi: twiddle factor angle.
   * @pxxi: first coefficient data memory address.
   * @pxxik: second coefficient data memory address.
   */
  int i, j, k, m, n, ncpy, step;
  real phi, *pxxi, *pxxik;

  /* extract the varargs. */
  int d = va_arg(*vl, int);
  real dir = (real) va_arg(*vl, double);
  hx_scalar *w = va_arg(*vl, hx_scalar*);
  hx_scalar *swp = va_arg(*vl, hx_scalar*);

  /* compute the number of bytes per scalar and the number of scalars. */
  ncpy = y->n * sizeof(real);
  n = y->sz[0];

  /* presort the scalar elements of the array. */
  for (i = 0, j = 0; i < n - 1; i++) {
    /* swap values according to sorting rules. */
    if (j > i) {
      /* swap all coefficients: x[i] <-> x[j]. */
      memcpy(swp->x, y->x + y->n * i, ncpy);
      memcpy(y->x + y->n * i, y->x + y->n * j, ncpy);
      memcpy(y->x + y->n * j, swp->x, ncpy);
    }

    /* right-shift the point count. */
    m = n >> 1;

    /* loop until this half of the array is traversed. */
    while (m <= j) {
      /* manipulation... */
      j -= m;
      m >>= 1;
    }

    /* increment the second-half array counter. */
    j += m;
  }

  /* initialize the transform outer loop counter. */
  k = 1;

  /* loop through the sorted data points. */
  do {
    /* compute the step value. */
    step = 2 * k;

    /* loop through the current segment of the array. */
    for (m = 0; m < k; m++) {
      /* compute the twiddle factor. */
      phi = M_PI * dir * (real) m / (real) k;
      hx_scalar_phasor(w, d, phi);

      /* loop through the other segment of the array. */
      for (i = m; i < n; i += step) {
        /* identify the memory addresses of the coefficients at
         * the (i) and (i+k) indices.
         */
        pxxi = y->x + y->n * i;
        pxxik = y->x + y->n * (i + k);

        /* compute the new values at the current memory locations:
         * swp <- w * x[i+k]
         * x[i+k] <- x[i] - swp;
         * x[i] <- x[i] + swp;
         */
        memset(swp->x, 0, ncpy);
        hx_data_mul(w->x, pxxik, swp->x, y->d, y->n, y->tbl);
        hx_data_add(pxxi, swp->x, pxxik, -1.0, y->d, y->n);
        hx_data_add(pxxi, swp->x, pxxi, 1.0, y->d, y->n);
      }
    }

    /* set (double) the loop counter value. */
    k = step;
  } while (k < n);

  /* return success. */
  return 1;
}

/* hx_array_fftfn(): computes an in-place radix-2 fast fourier transform along
 * a single dimension @d and direction @k of a hypercomplex multidimensional
 * array.
 * @x: pointer to the array structure.
 * @d: dimension to transform.
 * @k: direction to apply transform.
 * @dir: direction of transformation.
 */
int hx_array_fftfn (hx_array *x, int d, int k, real dir) {
  /* declare a few required variables:
   * @w: temporary array of twiddle factors.
   * @swp: temporary array of intermediate results, swapped values.
   * @nk: size of the transformation dimension.
   * @arr: index array used during iteration.
   * @xv: temporary vector for transformations.
   */
  hx_scalar w, swp;

  /* check that the dimensions are in bounds. */
  if (d < 0 || d >= x->d)
    throw("algebraic dimension %d out of bounds [0,%d)", d, x->d);

  /* check that the dimensions are in bounds. */
  if (k < 0 || k >= x->k)
    throw("topological dimension %d out of bounds [0,%d)", k, x->k);

  /* allocate temporary scalars for use in every transformation. */
  if (!hx_scalar_alloc(&w, x->d) ||
      !hx_scalar_alloc(&swp, x->d))
    throw("failed to allocate temporary %d-scalars", x->d);

  /* run the fft1d() callback function over every vector along @k */
  if (!hx_array_vector_op(x, k, &hx_array_fft1d, d, dir, &w, &swp))
    throw("failed to execute fft1d");

  /* free the temporary scalars. */
  hx_scalar_free(&w);
  hx_scalar_free(&swp);

  /* return success. */
  return 1;
}

/* hx_array_ht(): computes an in-place hilbert transform using a fast Fourier
 * transform to reconstruct the imaginary component of a signal from the
 * real component.
 * @x: pointer to the array structure.
 * @d: dimension to transform.
 * @k: direction to apply transform.
 */
int hx_array_ht (hx_array *x, int d, int k) {
  /* FIXME: implement hx_array_ht() */

  /* return success. */
  return 1;
}

/* hx_array_fshift(): circularly shifts each vector along a given array
 * topological dimension, but by a fractional amount of points.
 * @x: pointer to the array to manipulate.
 * @k: shift topological dimension.
 * @amount: fractional shift amount.
 */
int hx_array_fshift (hx_array *x, int k, real amount) {
  /* FIXME: implement hx_array_fshift() */

  /* return success. */
  return 1;
}

