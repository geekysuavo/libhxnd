
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

/* hx_array_ist_thresh(): soft-thresholds the values of one array into an
 * output array.
 * @xsrc: pointer to the source hypercomplex array.
 * @xdest: pointer to the destination hypercomplex array.
 * @lambda: pointer to the current thresholding magnitude, which should be
 *          zero on the first iteration to force an initialization.
 */
int hx_array_ist_thresh (hx_array *xsrc, hx_array *xdest, real *lambda) {
  /* declare a few required variables:
   * @i: hypercomplex scalar first coefficient index.
   * @xd: number of algebraic dimensions per scalar.
   * @xn: number of real coefficients per scalar.
   * @norm: current scalar norm value.
   */
  int i, xd, xn;
  real norm;

  /* locally store some array features. */
  xd = xsrc->d;
  xn = xsrc->n;

  /* check if the thresholding magnitude has been initialized. */
  if (*lambda <= 0.0) {
    /* loop over the elements to compute the maximum norm. */
    for (i = 0; i < xsrc->len; i += xn) {
      /* compute the current norm. */
      norm = hx_data_real_norm(xsrc->x + i, xn);

      /* check if the current norm exceeds the current maximum. */
      if (norm > *lambda)
        *lambda = norm;
    }
  }

  /* loop over the elements to apply the soft thresholding operation. */
  for (i = 0; i < xsrc->len; i += xn) {
    /* compute the current norm. */
    norm = hx_data_real_norm(xsrc->x + i, xn);

    /* determine whether the current norm exceeds the threshold. */
    if (norm > *lambda) {
      /* sum the element into the destination array. */
      hx_data_add(xdest->x + i, xsrc->x + i, xdest->x + i, 1.0, xd, xn);
    }
    else {
      /* zero the current source array element. */
      hx_data_zero(xsrc->x + i, xsrc->n);
    }
  }

  /* return success. */
  return 1;
}

/* hx_array_ist1d(): performs a set of one-dimensional reconstructions over a
 * two-dimensional array using iterative soft thresholding.
 * see hx_array_ist() for details.
 */
int hx_array_ist1d (hx_array *x, hx_index dx, hx_index kx,
                    int nsched, hx_index sched,
                    int niter, real thresh) {
  /* declare a few required variables:
   * @d: slice algebraic dimension index.
   * @k: slice topological dimension index.
   * @sz: current slice topological size.
   * @ja, @jb, @jmax: skipped iteration control variables.
   * @nzeros: number of unscheduled elements in @y.
   * @nbytes: number of bytes per hypercomplex scalar.
   * @zeros: linear indices of all unscheduled elements in @y.
   */
  int d, k, sz, ja, jb, jmax;
  int nzeros, nbytes;
  hx_index zeros;

  /* get the slice dimensionalities. */
  d = x->d;
  k = 1;

  /* get the slice length. */
  sz = x->sz[kx[1]];

  /* get the byte count per hypercomplex scalar. */
  nbytes = x->n * sizeof(real);

  /* get the number of norms and zeros. */
  nzeros = sz - nsched;

  /* allocate an array of unscheduled elements in each trace. */
  zeros = hx_index_unscheduled(k, &sz, 1, nsched, sched);

  /* check that allocation succeeded. */
  if (!zeros)
    throw("failed to allocate %d indices", nzeros);

  /* initialize the skipped iteration control variables. */
  hx_index_jump_init(x->k, x->sz, 1, &ja, &jb, &jmax);

  /* create a team of threads to execute multiple parallel reconstructions. */
  #pragma omp parallel
  {
    /* declare a few required thread-local variables:
     * @w: temporary twiddle-factor scalar.
     * @swp: temporary swap value scalar.
     * @y: currently sliced sub-array.
     * @z: output result sub-array.
     * @j: array skipped iteration master index.
     * @l: unscheduled array loop index.
     * @pidx: packed linear array index.
     * @iiter: main ist algorithm iteration loop counter.
     * @lambda: current iteration thresholding magnitude.
     */
    int j, l, pidx, iiter;
    hx_array y, yth, z;
    hx_scalar w, swp;
    real lambda;

    /* allocate temporary scalars for use in the fft. */
    if (!hx_scalar_alloc(&w, d) ||
        !hx_scalar_alloc(&swp, d))
      raise("failed to allocate temporary %d-scalars", d);

    /* allocate the slice destination arrays. */
    if (!hx_array_alloc(&z, d, k, &sz) ||
        !hx_array_alloc(&y, d, k, &sz) ||
        !hx_array_alloc(&yth, d, k, &sz))
      raise("failed to allocate temporary (%d, 1)-arrays", d);

    /* distribute tasks to the team of threads. */
    #pragma omp for
    for (j = 0; j < jmax; j++) {
      /* initialize the thresholding magnitude. */
      lambda = 0.0;

      /* compute the linear array index of the current vector. */
      pidx = hx_index_jump(j, ja, jb);

      /* slice the currently indexed vector from the array. */
      if (!hx_array_slice_vector(x, &y, 1, pidx))
        raise("failed to slice vector %d", j);

      /* loop over the iterations. */
      for (iiter = 0; iiter < niter; iiter++) {
        /* copy the current residual vector array. */
        memcpy(yth.x, y.x, y.len * sizeof(real));

        /* fourier transform the sliced vector array. */
        if (!hx_array_fft1d(&yth, dx[1], HX_FFT_FORWARD, &w, &swp))
          raise("failed to execute forward fft");

        /* threshold the frequency-domain vector. */
        if (!hx_array_ist_thresh(&yth, &z, &lambda))
          raise("failed to threshold sub-array %d", j);

        /* inverse fourier transform the sliced vector array. */
        if (!hx_array_fft1d(&yth, dx[1], HX_FFT_REVERSE, &w, &swp))
          raise("failed to execute inverse fft");

        /* reset the time-domain non-sampled points. */
        for (l = 0; l < nzeros; l++)
          memset(yth.x + yth.n * zeros[l], 0, nbytes);

        /* subtract the thresholded time-domain data from the residuals. */
        if (!hx_array_add_array(&y, &yth, -1.0, &y))
          raise("failed to deflate residual sub-array %d", j);

        /* scale down the threshold magnitude. */
        lambda *= thresh;
      }

      /* inverse fourier transform the reconstructed vector. */
      if (!hx_array_fft1d(&z, dx[1], HX_FFT_REVERSE, &w, &swp))
        raise("failed to execute final inverse fft");

      /* store the reconstructed vector back into the array. */
      if (!hx_array_store_vector(x, &z, 1, pidx))
        raise("failed to store vector %d", j);

      /* re-initialize the contents of the reconstructed vector. */
      hx_array_zero(&z);
    }

    /* free the temporary scalars. */
    hx_scalar_free(&w);
    hx_scalar_free(&swp);

    /* free the slice destination arrays. */
    hx_array_free(&z);
    hx_array_free(&y);
    hx_array_free(&yth);
  }

  /* free the array of unscheduled indices. */
  hx_index_free(zeros);

  /* return success. */
  return 1;
}

/* hx_array_istnd(): performs a set of multidimensional reconstructions over
 * an array having at least three dimensions using iterative soft thresholding.
 * see hx_array_ist() for details.
 */
int hx_array_istnd (hx_array *x, hx_index dx, hx_index kx,
                    int dsched, int nsched, hx_index sched,
                    int niter, real thresh) {
  /* declare a few required variables:
   * @iiter: iteration loop counter.
   * @lower: slice lower-bound index array.
   * @upper: slice upper-bound index array.
   * @ymax: maximum value in @ynorm.
   */
  hx_index sz, lower, upper;
  int i, j, n, d, k, iiter;

  /* @y: currently sliced sub-array.
   * @z: output result sub-array.
   */
  hx_array y, yth, z;

  /* @nzeros: number of unscheduled elements @y.
   * @zeros: linear indices of all unschedules elements in @y.
   * @lambda: current iteration thresholding magnitude.
   */
  int nbytes, nzeros;
  hx_index zeros;
  real lambda;

  /* get the number of reconstructions required. */
  n = x->sz[kx[0]];

  /* get the slice dimensionalities. */
  d = x->d;
  k = x->k;

  /* get the byte count per hypercomplex scalar. */
  nbytes = x->n * sizeof(real);

  /* allocate the slice size array. */
  sz = hx_index_alloc(k);

  /* check that allocation was successful. */
  if (!sz)
    throw("failed to allocate %d indices", k);

  /* build the slice size array. */
  for (i = 1, sz[0] = 1; i < k; i++)
    sz[i] = x->sz[kx[i]];

  /* allocate the bounding array indices. */
  lower = hx_index_alloc(k);
  upper = hx_index_alloc(k);

  /* check that the bounding arrays were allocated. */
  if (!lower || !upper)
    throw("failed to allocate bounding arrays");

  /* allocate the slice destination arrays. */
  if (!hx_array_alloc(&z, d, k, sz) ||
      !hx_array_alloc(&y, d, k, sz) ||
      !hx_array_alloc(&yth, d, k, sz))
    throw("failed to allocate (%d, %d)-arrays", d, k);

  /* determine the number of un-scheduled elements in each slice. */
  nzeros = y.len / y.n - nsched;
  zeros = hx_index_unscheduled(k - 1, sz + 1, dsched, nsched, sched);

  /* check that the array of zero indices was allocated. */
  if (!zeros)
    throw("failed to allocate %d indices", nzeros);

  /* store the elements of the bounding arrays. */
  for (i = 1; i < k; i++) {
    /* store the upper and lower bound. */
    upper[i] = x->sz[kx[i]] - 1;
    lower[i] = 0;
  }

  /* loop serially over the slices. */
  for (i = 0; i < n; i++) {
    /* initialize the thresholding magnitude. */
    lambda = 0.0;

    /* store the direct dimension bounds. */
    upper[0] = i;
    lower[0] = i;

    /* slice the indirect dimensions from the input array. */
    if (!hx_array_slice(x, &y, lower, upper))
      throw("failed to slice out sub-array %d", i);

    /* loop over the iterations. */
    for (iiter = 0; iiter < niter; iiter++) {
      /* copy the current residual data. */
      memcpy(yth.x, y.x, y.len * sizeof(real));

      /* loop over the sliced dimensions. */
      for (j = 1; j < k; j++) {
        /* forward fourier transform the slice. */
        if (!hx_array_fft(&yth, dx[j], kx[j]))
          throw("failed to apply forward fft");
      }

      /* threshold the frequency-domain sub-array. */
      if (!hx_array_ist_thresh(&yth, &z, &lambda))
        throw("failed to threshold sub-array %d", i);

      /* loop over the sliced dimensions. */
      for (j = 1; j < k; j++) {
        /* inverse fourier transform the slice. */
        if (!hx_array_ifft(&yth, dx[j], kx[j]))
          throw("failed to apply inverse fft");
      }

      /* reset the time-domain non-sampled points. */
      for (j = 0; j < nzeros; j++)
        memset(yth.x + yth.n * zeros[j], 0, nbytes);

      /* subtract the thresholded time-domain data from the residuals. */
      if (!hx_array_add_array(&y, &yth, -1.0, &y))
        throw("failed to deflate residual sub-array %d", i);

      /* scale down the threshold magnitude. */
      lambda *= thresh;
    }

    /* loop over the sliced dimensions. */
    for (j = 1; j < k; j++) {
      /* inverse fourier transform the reconstructed slice. */
      if (!hx_array_ifft(&z, dx[j], kx[j]))
        throw("failed to apply final inverse fft");
    }

    /* store the reconstructed slice back into the input array. */
    if (!hx_array_store(x, &z, lower, upper))
      throw("failed to store in sub-array %d", i);

    /* re-initialize the reconstructed slice for the next round. */
    hx_array_zero(&z);
  }

  /* free the slice destination arrays. */
  hx_array_free(&z);
  hx_array_free(&y);
  hx_array_free(&yth);

  /* free the allocated multidimensional indices. */
  hx_index_free(zeros);
  hx_index_free(lower);
  hx_index_free(upper);
  hx_index_free(sz);

  /* return success. */
  return 1;
}

/* hx_array_ist(): perform iterative soft thresholding to reconstruct
 * nonuniformly subsampled time-domain data in a hypercomplex array.
 * @x: pointer to the array to reconstruct.
 * @dx: array of algebraic dimension indices in @x.
 * @kx: array of topological dimension indices in @x.
 * @dsched: number of reconstruction dimensions.
 * @nsched: number of sampled @dsched-dimensional points.
 * @sched: array of @nsched @nnus-dimensional indices.
 * @niter: number of iterations to perform.
 * @thresh: threshold magnitude.
 */
int hx_array_ist (hx_array *x, hx_index dx, hx_index kx,
                  int dsched, int nsched, hx_index sched,
                  int niter, real thresh) {
  /* ensure the iteration count is in bounds. */
  if (niter < 1)
    throw("iteration count %d out of bounds [1,inf)", niter);

  /* ensure the threshold is in bounds. */
  if (thresh <= 0.0 || thresh >= 1.0)
    throw("threshold %.2f out of bounds (0,1)", thresh);

  /* determine which reconstruction function to use. */
  if (dsched == 1) {
    /* execute the one-dimensional function. */
    if (!hx_array_ist1d(x, dx, kx, nsched, sched, niter, thresh))
      throw("failed to execute one-dimensional reconstruction");
  }
  else if (dsched > 1) {
    /* execute the multidimensional function. */
    if (!hx_array_istnd(x, dx, kx, dsched, nsched, sched, niter, thresh))
      throw("failed to execute n-dimensional reconstruction");
  }
  else {
    /* throw an exception. */
    throw("invalid schedule configuration (%dx%d)", dsched, nsched);
  }

  /* return success. */
  return 1;
}

