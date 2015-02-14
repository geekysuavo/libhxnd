
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
 * of a hypercomplex one-dimensional (vector) array.
 * @y: pointer to the array structure to transform.
 * @d: dimension to transform.
 * @dir: direction of transformation.
 * @w: preallocated hypercomplex scalar.
 * @swp: preallocated hypercomplex scalar.
 */
int hx_array_fft1d (hx_array *y, int d, real dir,
                    hx_scalar *w, hx_scalar *swp) {
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
      phi = -M_PI * dir * (real) m / (real) k;
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

  /* for inverse transforms, scale each value in the vector. */
  if (dir == HX_FFT_REVERSE && !hx_array_scale(y, 1.0 / ((real) n), y))
    return 0;

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
   * @ja: small array stride for skipped iteration.
   * @jb: large array stride for skipped iteration.
   * @jmax: maximum loop control value.
   */
  int ja, jb, jmax;

  /* check that the dimensions are in bounds. */
  if (d < 0 || d >= x->d)
    throw("algebraic dimension %d out of bounds [0,%d)", d, x->d);

  /* check that the dimensions are in bounds. */
  if (k < 0 || k >= x->k)
    throw("topological dimension %d out of bounds [0,%d)", k, x->k);

  /* check that the transformation dimension is a power-of-two size. */
  if (!hx_ispow2(x->sz[k]))
    throw("dimension %d is not a power of two size (%d)", k, x->sz[k]);

  /* initialize the skipped iteration control variables. */
  hx_array_index_jump_init(x->k, x->sz, k, &ja, &jb, &jmax);

  /* create a team of threads to execute multiple parallel transforms. */
  #pragma omp parallel
  {
    /* declare a few required thread-local variables:
     * @j: array skipped iteration master index.
     * @idx: packed linear array index.
     * @w: temporary array of twiddle factors.
     * @swp: temporary array of intermediate results, swapped values.
     */
    int j, idx;
    hx_array xv;
    hx_scalar w, swp;

    /* allocate temporary scalars for use in every transformation. */
    if (!hx_scalar_alloc(&w, x->d) ||
        !hx_scalar_alloc(&swp, x->d))
      raise("failed to allocate temporary %d-scalars", x->d);

    /* allocate the slice destination vector array. */
    if (!hx_array_alloc(&xv, x->d, 1, &(x->sz[k])))
      raise("failed to allocate temporary (%d,1)-array", x->d);

    /* distribute tasks to the team of threads. */
    #pragma omp for
    for (j = 0; j < jmax; j++) {
      /* compute the linear array index of the current vector. */
      idx = hx_array_index_jump(j, ja, jb);

      /* slice the currently indexed vector from the array. */
      if (!hx_array_slice_vector(x, &xv, k, idx))
        raise("failed to slice vector %d", j);

      /* fourier transform the sliced vector array. */
      if (!hx_array_fft1d(&xv, d, dir, &w, &swp))
        raise("failed to execute vector fft %d", j);

      /* store the modified sliced vector back into the array. */
      if (!hx_array_store_vector(x, &xv, k, idx))
        raise("failed to store vector %d", j);
    }

    /* free the slice destination array. */
    hx_array_free(&xv);

    /* free the temporary scalars. */
    hx_scalar_free(&w);
    hx_scalar_free(&swp);
  }

  /* return success. */
  return 1;
}

/* hx_array_ht_cb(): callback function for hx_array_ht().
 *
 * args:
 *  see hx_array_foreach_cb().
 *
 * varargs:
 *  @xtmp: temporary swap vector.
 */
int hx_array_ht_cb (hx_array *x, hx_array *y,
                    int *arr, int idx,
                    va_list *vl) {
  /* declare required variables. */
  int n, off0, off1;

  /* extract the varargs. */
  hx_array *xtmp = va_arg(*vl, hx_array*);

  /* get the size of the vector. */
  n = y->sz[0];

  /* copy the slice into the temporary location. */
  memcpy(xtmp->x, y->x, y->len * sizeof(real));

  /* compute offsets to the (i=0) and (i=n/2) coefficients. */
  off0 = 0;
  off1 = (n / 2) * y->n;

  /* build up the even-length shuffled result:
   *  [y(0), 2 * y(1 : n/2 - 1), y(n/2), zeros(n/2-1)]
   *   0     1 .. n/2 - 1        n/2
   */
  hx_array_zero(y);
  hx_data_add(NULL, xtmp->x + off0, y->x + off0, 0.5, y->d, y->n);
  hx_data_add(NULL, xtmp->x + off1, y->x + off1, 0.5, y->d, y->n);
  memcpy(y->x + y->n, xtmp->x + y->n, ((n / 2) - 1) * y->n * sizeof(real));

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
  /* declare a few required variables:
   * @xtmp: temporary duplicate array of @x.
   */
  hx_array xtmp;
  int i, n;

  /* check that the algebraic dimension index is in bounds. */
  if (d < 0 || d >= x->d)
    throw("transform index %d out of bounds [0,%d)", d, x->d);

  /* check that the topological dimension index is in bounds. */
  if (k < 0 || k >= x->k)
    throw("shift index %d out of bounds [0,%d)", k, x->k);

  /* get the size of the transformation dimension. */
  n = x->sz[k];

  /* allocate a temporary duplicate array. */
  if (!hx_array_alloc(&xtmp, x->d, 1, &n))
    throw("failed to allocate temporary array");

  /* zero the d-dimension values in the array. */
  for (i = 0; i < xtmp.len; i += xtmp.n)
    xtmp.x[(1 << d) + i] = 0.0;

  /* forward Fourier-transform the vectors. */
  if (!hx_array_fft(x, d, k))
    throw("failed to apply forward fft");

  /* perform the data shuffling. */
  if (!hx_array_foreach_vector(x, k, &hx_array_ht_cb, &xtmp))
    throw("failed to apply shuffling operation");

  /* inverse Fourier-transform the vectors. */
  if (!hx_array_ifft(x, d, k))
    throw("failed to apply inverse fft");

  /* free the temporary duplicate array. */
  hx_array_free(&xtmp);

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
  /* declare a few required variables:
   * @phi: scalar phasor entries in @ph.
   * @ph: linear phase adjustment array.
   * @i: general-purpose loop counter.
   * @fi: fractional loop counter.
   * @n: shift dimension size.
   */
  hx_scalar phi;
  hx_array ph;
  int i, j, n;
  real fi;

  /* check that the shift dimension index is in bounds. */
  if (k < 0 || k >= x->k)
    throw("shift index %d out of bounds [0,%d)", k, x->k);

  /* check if the shift amount is zero. */
  if (amount == 0.0)
    return 1;

  /* locally store the size of the shift dimension. */
  n = x->sz[k];

  /* allocate the linear phase array. */
  if (!hx_scalar_alloc(&phi, x->d) ||
      !hx_array_alloc(&ph, x->d, 1, &n))
    throw("failed to allocate temporary phase array");

  /* compute the values that will reside in the linear phase array. */
  for (i = 0; i < n; i++) {
    /* compute the output index. */
    if (n % 2)
      j = (i > n / 2 ? i - n / 2 - 1 : i + n / 2);
    else
      j = (i < n / 2 ? i + n / 2 : i - n / 2);

    /* compute the fractional loop index. */
    fi = 2.0 * ((real) i) / ((real) (n - 1)) - 1.0;

    /* compute the scalar phasing element. */
    hx_scalar_phasor(&phi, k, -M_PI * amount * fi);

    /* copy the phasor scalar element into the array. */
    memcpy(ph.x + ph.n * j, phi.x, ph.n * sizeof(real));
  }

  /* forward Fourier-transform the vectors. */
  if (!hx_array_fft(x, k, k))
    throw("failed to apply forward fft");

  /* perform the scaling. */
  if (!hx_array_mul_vector(x, &ph, k, x))
    throw("failed to apply linear phase");

  /* inverse Fourier-transform the vectors. */
  if (!hx_array_ifft(x, k, k))
    throw("failed to apply inverse fft");

  /* free the linear phase array. */
  hx_scalar_free(&phi);
  hx_array_free(&ph);

  /* return success. */
  return 1;
}

