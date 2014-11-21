
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
 * @norms: pre-allocated array to hold array norm values.
 * @thresh: relative threshold value to apply.
 */
int hx_array_ist_thresh (hx_array *xsrc, hx_array *xdest,
                         real *norms, real thresh) {
  /* declare a few required variables. */
  int i, j, n, xd, xn;
  real maxnorm, tau;

  /* determine the number of norm elements. */
  n = xsrc->len / xsrc->n;

  /* locally store some array features. */
  xd = xsrc->d;
  xn = xsrc->n;

  /* compute the 'inverse' of the threshold value. */
  tau = -(1.0 - thresh);

  /* loop over the elements to compute the norms. */
  for (i = 0, maxnorm = 0.0; i < n; i++) {
    /* compute the norm. */
    norms[i] = hx_data_real_norm(xsrc->x + i * xsrc->n, xsrc->d, xsrc->n);

    /* see if the new norm is the largest so far. */
    if (norms[i] > maxnorm)
      maxnorm = norms[i];
  }

  /* scale the maximum norm value by the threshold. */
  maxnorm *= thresh;

  /* loop a final time over the elements to apply the threshold. */
  for (i = 0, j = 0; i < n; i++, j += xn) {
    /* determine whether the current norm exceeds the threshold. */
    if (norms[i] > maxnorm) {
      /* sum the element into the destination array. */
      hx_data_add(xdest->x + j, xsrc->x + j, xdest->x + j, 1.0, xd, xn);

      /* subtract the thresholded element from the source array. */
      hx_data_add(xsrc->x + j, xsrc->x + j, xsrc->x + j, tau, xd, xn);
    }
  }

  /* return success. */
  return 1;
}

/* hx_array_ist1d(): performs a set of one-dimensional reconstructions over a
 * two-dimensional array using iterative soft thresholding.
 * see hx_array_ist() for details.
 */
int hx_array_ist1d (hx_array *x, int nsched, int *sched,
                    int niter, real thresh) {
  /* declare a few required variables:
   * @iiter: iteration loop counter.
   * @idx: slice location.
   */
  int iiter, ja, jb, jmax;

  /* @nnorms: number of scalar elements in @y.
   * @nzeros: number of unscheduled elements in @y.
   * @nbytes: number of bytes per hypercomplex scalar.
   */
  int d, k, sz, nnorms, nzeros, nbytes, *zeros;

  /* get the slice dimensionalities. */
  d = x->d;
  k = 1;

  /* get the slice length. */
  sz = x->sz[1];

  /* get the byte count per hypercomplex scalar. */
  nbytes = x->n * sizeof(real);

  /* get the number of norms and zeros. */
  nnorms = sz;
  nzeros = sz - nsched;

  /* allocate an array of unscheduled elements in each trace. */
  zeros = hx_array_index_unscheduled(k, &sz, 1, nsched, sched);

  /* check that allocation succeeded. */
  if (!zeros)
    throw("failed to allocate %d indices", nzeros);

  /* initialize the skipped iteration control variables. */
  hx_array_index_jump_init(x->k, x->sz, 1, &ja, &jb, &jmax);

  /* create a team of threads to execute multiple parallel reconstructions. */
  #pragma omp parallel
  {
    /* declare a few required thread-local variables:
     * @tid: thread index.
     * @tnum: thread count.
     */
    int j, l, idx, tid, tnum;

    /* @y: currently sliced sub-array.
     * @z: output result sub-array.
     */
    hx_scalar w, swp;
    hx_array y, z;

    /* @ynorms: norms of all values in @y.
     * @zeros: linear indices of all unscheduled elements in @y.
     */
    real *ynorms;

    /* gain a handle of the thread index and number of threads. */
    tnum = omp_get_num_threads();
    tid = omp_get_thread_num();

    /* allocate temporary scalars for use in the fft. */
    if (!hx_scalar_alloc(&w, d) ||
        !hx_scalar_alloc(&swp, d))
      raise("failed to allocate temporary %d-scalars", d);

    /* allocate the slice destination arrays. */
    if (!hx_array_alloc(&y, d, k, &sz) ||
        !hx_array_alloc(&z, d, k, &sz))
      raise("failed to allocate temporary (%d,1)-arrays", d);

    /* allocate an array of d-norm values. */
    ynorms = (real*) calloc(nnorms, sizeof(real));

    /* distribute tasks to the team of threads. */
    for (j = tid; j < jmax; j += tnum) {
      /* compute the linear array index of the current vector. */
      idx = hx_array_index_jump(j, ja, jb);

      /* slice the currently indexed vector from the array. */
      if (!hx_array_slice_vector(x, &y, 1, idx))
        raise("failed to slice vector %d", j);

      /* loop over the iterations. */
      for (iiter = 0; iiter < niter; iiter++) {
        /* fourier transform the sliced vector array. */
        if (!hx_array_fft1d(&y, 1, HX_FFT_FORWARD, &w, &swp))
          raise("failed to execute forward fft");

        /* threshold the frequency-domain vector. */
        if (!hx_array_ist_thresh(&y, &z, ynorms, thresh))
          raise("failed to threshold sub-array %d", j);

        /* inverse fourier transform the sliced vector array. */
        if (!hx_array_fft1d(&y, 1, HX_FFT_REVERSE, &w, &swp))
          raise("failed to execute inverse fft");

        /* reset the time-domain non-sampled points. */
        for (l = 0; l < nzeros; l++)
          memset(y.x + y.n * zeros[l], 0, nbytes);
      }

      /* store the modified sliced vector back into the array. */
      if (!hx_array_store_vector(x, &z, 1, idx))
        raise("failed to store vector %d", j);
    }

    /* free the d-norm array. */
    free(ynorms);

    /* free the temporary scalars. */
    hx_scalar_free(&w);
    hx_scalar_free(&swp);

    /* free the slice destination arrays. */
    hx_array_free(&y);
    hx_array_free(&z);
  }

  /* free the array of unscheduled indices. */
  free(zeros);

  /* return success. */
  return 1;
}

/* hx_array_istnd(): performs a set of multidimensional reconstructions over
 * an array having at least three dimensions using iterative soft thresholding.
 * see hx_array_ist() for details.
 */
int hx_array_istnd (hx_array *x, int dsched, int nsched, int *sched,
                    int niter, real thresh) {
  /* declare a few required variables:
   * @iiter: iteration loop counter.
   * @lower: slice lower-bound index array.
   * @upper: slice upper-bound index array.
   * @ymax: maximum value in @ynorm.
   */
  int i, j, n, d, k, *sz;
  int iiter, *lower, *upper;

  /* @y: currently sliced sub-array.
   * @z: output result sub-array.
   */
  hx_array y, z;

  /* @nnorms: number of scalar elements in @y.
   * @ynorms: norms of all values in @y.
   */
  real *ynorms;
  int nnorms;

  /* @nzeros: number of unscheduled elements @y.
   * @zeros: linear indices of all unschedules elements in @y.
   */
  int nbytes, nzeros, *zeros;

  /* get the number of reconstructions required. */
  n = x->sz[0];

  /* get the slice dimensionalities. */
  d = x->d;
  k = x->k;

  /* get the byte count per hypercomplex scalar. */
  nbytes = x->n * sizeof(real);

  /* allocate the slice size array. */
  sz = hx_array_index_alloc(k);

  /* check that allocation was successful. */
  if (!sz)
    throw("failed to allocate %d indices", k);

  /* build the slice size array. */
  for (i = 1, sz[0] = 1; i < k; i++)
    sz[i] = x->sz[i];

  /* allocate the bounding array indices. */
  lower = hx_array_index_alloc(k);
  upper = hx_array_index_alloc(k);

  /* check that the bounding arrays were allocated. */
  if (!lower || !upper)
    throw("failed to allocate bounding arrays");

  /* allocate the slice destination arrays. */
  if (!hx_array_alloc(&y, d, k, sz) ||
      !hx_array_alloc(&z, d, k, sz))
    throw("failed to allocate (%d, %d)-array", d, k);

  /* allocate an array of d-norm values. */
  nnorms = y.len / y.n;
  ynorms = (real*) calloc(nnorms, sizeof(real));

  /* check that the array of d-norm values was allocated. */
  if (!ynorms)
    throw("failed to allocate array of %d norms", nnorms);

  /* determine the number of un-scheduled elements in each slice. */
  nzeros = nnorms - nsched;
  zeros = hx_array_index_unscheduled(k - 1, sz + 1, dsched, nsched, sched);

  /* check that the array of zero indices was allocated. */
  if (!zeros)
    throw("failed to allocate %d indices", nzeros);

  /* build the bounding arrays. */
  for (i = 1; i < k; i++) {
    /* store the upper and lower bound. */
    upper[i] = x->sz[i] - 1;
    lower[i] = 0;
  }

  /* loop serially over the slices. */
  for (i = 0; i < n; i++) {
    /* store the direct dimension bounds. */
    upper[0] = i;
    lower[0] = i;

    /* slice the indirect dimensions from the input array. */
    if (!hx_array_slice(x, &y, lower, upper))
      throw("failed to slice out sub-array %d", i);

    /* loop over the iterations. */
    for (iiter = 0; iiter < niter; iiter++) {
      /* loop over the sliced dimensions. */
      for (j = 1; j < k; j++) {
        /* forward Fourier-transform the slice. */
        if (!hx_array_fft(&y, j, j))
          throw("failed to apply forward fft");
      }

      /* threshold the frequency-domain sub-array. */
      if (!hx_array_ist_thresh(&y, &z, ynorms, thresh))
        throw("failed to threshold sub-array %d", i);

      /* loop over the sliced dimensions. */
      for (j = 1; j < k; j++) {
        /* inverse Fourier-transform the slice. */
        if (!hx_array_ifft(&y, j, j))
          throw("failed to apply inverse fft");
      }

      /* reset the time-domain non-sampled points. */
      for (j = 0; j < nzeros; j++)
        memset(y.x + y.n * zeros[j], 0, nbytes);
    }

    /* store the reconstructed slice back into the input array. */
    if (!hx_array_store(x, &z, lower, upper))
      throw("failed to store in sub-array %d", i);
  }

  /* free the array of un-scheduled indices. */
  free(zeros);

  /* free the array of norm values. */
  free(ynorms);

  /* free the slice destination arrays. */
  hx_array_free(&y);
  hx_array_free(&z);

  /* free the allocated index arrays. */
  free(lower);
  free(upper);
  free(sz);

  /* return success. */
  return 1;
}

/* hx_array_ist(): performs iterative soft thresholding to reconstruct
 * nonuniformly subsampled time-domain data in a hypercomplex array.
 * @x: pointer to the array to reconstruct.
 * @dsched: number of reconstruction dimensions.
 * @nsched: number of sampled @dsched-dimensional points.
 * @sched: array of @nsched @nnus-dimensional indices.
 * @niter: number of iterations to perform.
 * @thresh: threshold magnitude.
 */
int hx_array_ist (hx_array *x, int dsched, int nsched, int *sched,
                  int niter, real thresh) {
  /* determine which reconstruction function to use. */
  if (dsched == 1) {
    /* execute the one-dimensional function. */
    return hx_array_ist1d(x, nsched, sched, niter, thresh);
  }
  else if (dsched > 1) {
    /* execute the multidimensional function. */
    return hx_array_istnd(x, dsched, nsched, sched, niter, thresh);
  }

  /* return failure, because @dsched < 1 */
  throw("invalid schedule dimensionality");
}
