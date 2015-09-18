
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

/* include the n-dimensional math header. */
#include <hxnd/hx.h>

/* hx_array_ffm1d(): perform a set of one-dimensional fast forward
 * maximum entropy reconstructions over a two-dimensional array.
 * see hx_array_ffm() for details.
 */
int hx_array_ffm1d (hx_array *x, hx_index dx, hx_index kx,
                    int nsched, hx_index sched, int niter,
                    hx_entropy_functional f,
                    hx_entropy_functional df) {
  /* declare a few required variables:
   * @d: slice algebraic dimension index.
   * @k: slice topological dimension index.
   * @sz: twice the current slice topological size.
   * @ja, @jb, @jmax: skipped iteration control variables.
   * @nzeros: number of unscheduled elements in @xj.
   * @nbytes: number of bytes per hypercomplex scalar.
   * @zeros: linear indices of all unscheduled elements in @xj.
   * @alpha: fixed iteration step scale factor.
   */
  int d, k, sz, ja, jb, jmax;
  int nbytes;
  real alpha;

  /* get the slice dimensionalities. */
  d = x->d;
  k = 1;

  /* get the slice length. */
  sz = 2 * x->sz[kx[1]];
  alpha = (real) sz;

  /* get the byte count per hypercomplex scalar. */
  nbytes = x->n * sizeof(real);

  /* initialize the skipped iteration control variables. */
  hx_index_jump_init(x->k, x->sz, kx[1], &ja, &jb, &jmax);

  /* create a team of threads to execute multiple parallel reconstructions. */
  #pragma omp parallel
  {
    /* declare a few thread-local variables:
     * @j: array skipped iteration master index.
     * @l: unscheduled array loop index.
     * @pidx: packed linear array index.
     * @iiter: ffm iteration loop counter.
     * @xj: currently sliced sub-array.
     * @g: current gradient sub-array.
     * @w: temporary twiddle-factor scalar.
     * @swp: temporary swap value scalar.
     * @beta: conjugate gradient step factor.
     */
    int j, l, pidx, iiter;
    hx_scalar w, swp;
    hx_array xj, g;

    /* allocate temporary scalars to use in the fft. */
    if (!hx_scalar_alloc(&w, d) ||
        !hx_scalar_alloc(&swp, d))
      raise("failed to allocate temporary %d-scalars", d);

    /* allocate the scratch-space arrays. */
    if (!hx_array_alloc(&g, d, k, &sz) ||
        !hx_array_alloc(&xj, d, k, &sz))
      raise("failed to allocate temporary (%d, 1)-arrays", d);

    /* distribute tasks to the team of threads. */
    #pragma omp for
    for (j = 0; j < jmax; j++) {
      /* compute the linear array index of the current vector. */
      pidx = hx_index_jump(j, ja, jb);

      /* slice the currently indexed vector from the array. */
      if (!hx_array_slice_vector(x, &xj, kx[1], pidx))
        raise("failed to slice vector %d", j);

      /* loop over the iterations. */
      for (iiter = 0; iiter < niter; iiter++) {
        /* initialize the gradient vector. */
        memcpy(g.x, xj.x, xj.len * sizeof(real));

        /* fourier transform the time-domain vector. */
        if (!hx_array_fft1d(&g, dx[1], HX_FFT_FORWARD, &w, &swp))
          raise("failed to execute forward fft");

        /* compute the gradient of the entropy. */
        for (l = 0; l < g.len; l += g.n)
          df(g.x + l, g.x + l, g.n);

        /* inverse fourier transform the gradient vector. */
        if (!hx_array_fft1d(&g, dx[1], HX_FFT_REVERSE, &w, &swp))
          raise("failed to execute inverse fft");

        /* reset the sampled time-domain points in the gradient. */
        for (l = 0; l < nsched; l++)
          memset(g.x + g.n * sched[l], 0, nbytes);

        /* update the time-domain vector. */
        if (!hx_array_add_array(&xj, &g, alpha, &xj))
          raise("failed to update time-domain vector");
      }

      /* store the reconstructed vector back into the array. */
      if (!hx_array_store_vector(x, &xj, kx[1], pidx))
        raise("failed to store vector %d", j);
    }

    /* free the temporary scalars. */
    hx_scalar_free(&w);
    hx_scalar_free(&swp);

    /* free the scratch-space arrays. */
    hx_array_free(&g);
    hx_array_free(&xj);
  }

  /* return success. */
  return 1;
}

/* hx_array_ffmnd(): perform a set of multidimensional fast forward
 * maximum entropy reconstructions over an array having at least
 * three dimensions.
 * see hx_array_ffm() for details.
 */
int hx_array_ffmnd (hx_array *x, hx_index dx, hx_index kx,
                    int dsched, int nsched, hx_index sched, int niter,
                    hx_entropy_functional f,
                    hx_entropy_functional df) {
  /* declare a few required variables:
   * @sz: size of the temporary arrays.
   * @lower: slice lower-bound index array.
   * @upper: slice upper-bound index array.
   * @iiter: iteration loop counter.
   */
  hx_index sz, lower, upper;
  int i, j, n, d, k, iiter;

  /* @xi: currently sliced sub-array.
   * @g: current gradient sub-array.
   */
  hx_array xi, g;

  /* @nbytes: number of bytes per hypercomplex scalar.
   * @alpha: fixed iteration step scale factor.
   */
  real alpha;
  int nbytes;

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
  for (i = 0, sz[0] = 1; i < k; i++)
    sz[i] = 2 * x->sz[kx[i]];

  /* allocate the bounding array indices. */
  lower = hx_index_alloc(k);
  upper = hx_index_alloc(k);

  /* check that the bounding arrays were allocated. */
  if (!lower || !upper)
    throw("failed to allocate bounding arrays");

  /* allocate the scratch-space arrays. */
  if (!hx_array_alloc(&g, d, k, sz) ||
      !hx_array_alloc(&xi, d, k, sz))
    throw("failed to allocate temporary (%d, %d)-arrays", d, k);

  /* store the elements of the bounding arrays. */
  for (i = 1; i < k; i++) {
    /* store the upper and lower bound. */
    upper[i] = x->sz[kx[i]] - 1;
    lower[i] = 0;
  }

  /* compute the scale factor. */
  for (i = 1, alpha = 1.0; i < k; i++)
    alpha *= (real) sz[i];

  /* loop serially over the slices. */
  for (i = 0; i < n; i++) {
    /* store the direct dimension bounds. */
    upper[0] = i;
    lower[0] = i;

    /* slice the indirect dimensions from the input array. */
    if (!hx_array_slice(x, &xi, lower, upper))
      throw("failed to slice out sub-array %d", i);

    /* loop over the iterations. */
    for (iiter = 0; iiter < niter; iiter++) {
      /* initialize the gradient array. */
      memcpy(g.x, xi.x, xi.len * sizeof(real));

      /* loop over the sliced dimensions. */
      for (j = 1; j < k; j++) {
        /* forward fourier transform the slice. */
        if (!hx_array_fft(&g, dx[j], kx[j]))
          throw("failed to apply forward fft");
      }

      /* compute the gradient of the entropy. */
      for (j = 0; j < g.len; j += g.n)
        df(g.x + j, g.x + j, g.n);

      /* loop over the sliced dimensions. */
      for (j = 1; j < k; j++) {
        /* inverse fourier transform the slice. */
        if (!hx_array_ifft(&g, dx[j], kx[j]))
          throw("failed to apply forward fft");
      }

      /* reset the sampled time-domain points in the gradient */
      for (j = 0; j < nsched; j++)
        memset(g.x + g.n * sched[j], 0, nbytes);

      /* update the time-domain vector. */
      if (!hx_array_add_array(&xi, &g, alpha, &xi))
        throw("failed to update time-domain array");
    }

    /* store the reconstructed slice back into the input array. */
    if (!hx_array_store(x, &xi, lower, upper))
      throw("failed to store in sub-array %d", i);
  }

  /* free the scratch-space arrays. */
  hx_array_free(&g);
  hx_array_free(&xi);

  /* free the allocated multidimensional indices. */
  hx_index_free(lower);
  hx_index_free(upper);
  hx_index_free(sz);

  /* return success. */
  return 1;
}

/* hx_array_ffm(): use the fast forward maximum entropy algorithm to
 * reconstruct nonuniformly subsampled time-domain data in a hypercomplex
 * array.
 * @x: pointer to the array to reconstruct.
 * @dx: array of algebraic dimension indices in @x.
 * @kx: array of topological dimension indices in @x.
 * @dsched: number of reconstruction dimensions.
 * @nsched: number of sampled @dsched-dimensional points.
 * @niter: maximum number of iterations to perform.
 * @type: entropy functional type to utilize..
 */
int hx_array_ffm (hx_array *x, hx_index dx, hx_index kx,
                  int dsched, int nsched, hx_index sched,
                  int niter, enum hx_entropy_type type) {
  /* declare a few required variables:
   * @f: function pointer for the entropy functional.
   * @df: function pointer for the entropy derivative.
   */
  hx_entropy_functional f, df;

  /* ensure the schedule is allocated and properly sized. */
  if (!sched || dsched < 1 || nsched < 1)
    throw("invalid schedule configuration (%dx%d)", dsched, nsched);

  /* ensure the iteration count is in bounds. */
  if (niter < 1)
    throw("iteration count %d out of bounds [1,inf)", niter);

  /* retrieve the entropy functional function pointers. */
  if (!hx_entropy_get_functionals(type, &f, &df))
    throw("failed to retrieve entropy functionals");

  /* determine which reconstruction function to use. */
  if (dsched == 1) {
    /* execute the one-dimensional function. */
    if (!hx_array_ffm1d(x, dx, kx, nsched, sched, niter, f, df))
      throw("failed to execute one-dimensional reconstruction");
  }
  else if (dsched > 1) {
    /* execute the multidimensional function. */
    if (!hx_array_ffmnd(x, dx, kx, dsched, nsched, sched, niter, f, df))
      throw("failed to execute n-dimensional reconstruction");
  }
  else {
    /* throw an exception. */
    throw("invalid schedule configuration (%dx%d)", dsched, nsched);
  }

  /* return success. */
  return 1;
}

