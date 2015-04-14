
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

/* hx_array_irls_dftmatrix(): compute the matrix form of the multidimensional
 * inverse discrete Fourier transform (IDFT) matrix of a nonuniformly sampled
 * time-domain array.
 * @d: number of array algebraic dimensions.
 * @k: number of array topological dimensions.
 * @n: number of scheduled time-domain indices.
 * @sz: topological sizes of the array.
 * @dx: array of algebraic transform dimensions.
 * @kx: array of topological transform dimensions.
 * @sched: linear array of packed scheduled indices.
 * @F: pointer to the output hypercomplex array.
 */
int hx_array_irls_dftmatrix (int d, int k, int n, hx_index sz,
                             hx_index dx, hx_index kx,
                             hx_index sched,
                             hx_array *F) {
  /* declare required variables:
   * @Fsz: sizes of the output array.
   * @Fidx: indices inside the output array.
   * @idxi: multidimensional time-domain index.
   * @idxk: multidimensional frequency-domain index.
   * @ph: current output matrix element value.
   * @phd: current dimension phase factor.
   * @phtmp: temporary phase factor.
   * @i, @j, @di: general-purpose loop counters.
   * @N: number of output frequency-domain indices.
   * @Fpidx: packed linear output matrix index.
   * @theta: current dimension phase angle.
   */
  hx_index Fsz, Fidx, idxi, idxk;
  hx_scalar ph, phd, phtmp;
  int i, j, di, N, Fpidx;
  real theta;

  /* compute the number of frequency-domain indices. */
  for (i = 1, N = 1; i < k; i++)
    N *= sz[kx[i]];

  /* allocate the input and output multidimensional indices. */
  idxi = hx_index_alloc(k);
  idxk = hx_index_alloc(k);

  /* allocate the output array size index. */
  Fsz = hx_index_build(2, n, N);
  Fidx = hx_index_alloc(2);

  /* ensure the indices were allocated. */
  if (!idxi || !idxk || !Fsz || !Fidx)
    throw("failed to allocate multidimensional indices");

  /* allocate temporary hypercomplex scalars. */
  if (!hx_scalar_alloc(&ph, d) ||
      !hx_scalar_alloc(&phd, d) ||
      !hx_scalar_alloc(&phtmp, d))
    throw("failed to allocate %d-scalars", d);

  /* ensure the output array is allocated properly. */
  if (F->k != 2 || F->d != d || hx_index_cmp(F->k, F->sz, Fsz)) {
    /* nope. attempt to allocate it. */
    if (!hx_array_alloc(F, d, 2, Fsz))
      throw("failed to allocate output (%d, 2)-array", d);
  }

  /* loop over the output matrix rows. */
  for (i = 0; i < n; i++) {
    /* unpack the current scheduled time-domain index. */
    hx_index_unpack(k - 1, sz + 1, idxi, sched[i]);

    /* loop over the output matrix columns. */
    for (j = 0; j < N; j++) {
      /* unpack the current frequency-domain index. */
      hx_index_unpack(k - 1, sz + 1, idxk, j);

      /* initialize the final scalar phase factor. */
      hx_scalar_zero(&ph);
      ph.x[0] = 1.0;

      /* loop over the transform dimensions. */
      for (di = 1; di < k; di++) {
        /* compute the current dimension phase angle. */
        theta = 2.0 * M_PI;
        theta *= (real) idxi[kx[di]];
        theta *= (real) idxk[kx[di]];
        theta /= (real) sz[kx[di]];

        /* 1. compute the current dimension phase factor.
         * 2. copy the current phase factor into a temporary scalar.
         * 3. multiply the phase factor into the final phase factor.
         */
        if (!hx_data_copy(ph.x, phtmp.x, ph.n) ||
            !hx_scalar_phasor(&phd, dx[di], theta) ||
            !hx_scalar_mul(&phtmp, &phd, &ph))
          throw("failed to compute phase factors");
      }

      /* scale the phase factor into the final matrix element. */
      if (!hx_scalar_scale(&ph, 1.0 / sqrt((real) N), &ph))
        throw("failed to scale phase factor");

      /* compute the output matrix element packed linear index. */
      Fidx[0] = i;
      Fidx[1] = j;
      hx_index_pack(2, Fsz, Fidx, &Fpidx);

      /* store the output matrix element. */
      if (!hx_data_copy(ph.x, F->x + Fpidx * F->n, F->n))
        throw("failed to store phase factor");
    }
  }

  /* free the allocated hypercomplex scalars. */
  hx_scalar_free(&phtmp);
  hx_scalar_free(&phd);
  hx_scalar_free(&ph);

  /* free the allocated indices. */
  hx_index_free(idxi);
  hx_index_free(idxk);
  hx_index_free(Fidx);
  hx_index_free(Fsz);

  /* return success. */
  return 1;
}

/* hx_array_irlsfun(): compute a single irls reconstruction.
 * @X: vector of frequency-domain spectral estimates.
 * @x: vector of time-domain acquired values.
 * @w: vector of real spectral weights.
 * @z: vector of time-domain estimates.
 * @wx: vector of weighted spectral estimates.
 * @rx: vector of time-domain residuals.
 * @A: hermitian positive definite array of regression factors.
 * @niter, @pa, @pb: see hx_array_irls().
 */
int hx_array_irlsfun (hx_array *X, hx_array *x, hx_array *w, hx_array *z,
                      hx_array *wx, hx_array *rx, hx_array *A,
                      int niter, real pa, real pb) {
  /* declare a few required variables:
   */
  real p, dp;
  int iiter;

  /* compute the change in norm p-value per iteration. */
  dp = (pb - pa) / (real) (niter - 1);

  /* FIXME: compute the initial spectral estimate. */

  /* loop over the reconstruction iterations. */
  for (iiter = 0; iiter < niter; iiter++) {
    /* compute the current norm p-value. */
    p = pa + (real) iiter * dp;

    /* FIXME: compute the current weights. */

    /* FIXME: compute the current weighted spectral estimate. */

    /* FIXME: compute the current time-domain residual. */

    /* FIXME: compute the sum of squared weighted estimates. */

    /* FIXME: compute the sum of squared residuals. */

    /* FIXME: compute the current lagrange multiplier value. */

    /* FIXME: compute the final inverted weights. */

    /* FIXME: compute the normal matrix. */

    /* FIXME: decompose the normal matrix. */

    /* FIXME: solve for the time-domain factors. */

    /* FIXME: compute the new spectral estimate. */
  }

  /* return success. */
  return 1;
}

/* hx_array_irls(): perform iteratively reweighted least squares regression
 * to reconstruct nonuniformly subsampled time-domain data in a hypercomplex
 * array.
 * @x: pointer to the array to reconstruct.
 * @dx: array of algebraic dimension indices in @x.
 * @kx: array of topological dimension indices in @x.
 * @dsched: number of reconstruction dimensions.
 * @nsched: number of sampled @dsched-dimensional points.
 * @niter: number of iterations to perform.
 * @pa: starting norm p-value.
 * @pb: ending norm p-value.
 */
int hx_array_irls (hx_array *x, hx_index dx, hx_index kx,
                   int dsched, int nsched, hx_index sched,
                   int niter, real pa, real pb) {
  /* declare a few required variables:
   * @F: discrete Fourier transform matrix.
   * @Y: spectral estimate array.
   * @y: time-domain data vector.
   * @w: weight vector.
   * @psched: packed schedule indices.
   */
  hx_index Asz, sz, xsched, ysched, lower, upper;
  hx_array A, F, Y, y, w, z, wx, rx;
  int d, k, i, n, is, ns, N;

  /* ensure the schedule is allocated and properly sized. */
  if (!sched || dsched < 1 || nsched < 1)
    throw("invalid schedule configuration (%dx%d)", dsched, nsched);

  /* ensure the iteration count is in bounds. */
  if (niter < 1)
    throw("iteration count %d out of bounds [1,inf)", niter);

  /* ensure the starting p-value is in bounds. */
  if (pa < 0.0 || pa > 1.0)
    throw("starting norm order %.3d out of bounds [0,1]", pa);

  /* ensure the ending p-value is in bounds. */
  if (pb < 0.0 || pb > 1.0)
    throw("ending norm order %.3d out of bounds [0,1]", pb);

  /* ensure the p-values are decreasing. */
  if (pa < pb)
    throw("norm orders must decrease during iteration");

  /* store the input array dimensionalities. */
  d = x->d;
  k = x->k;

  /* store the number of slices to reconstruct. */
  ns = x->sz[kx[0]];

  /* allocate a size multidimensional index. */
  sz = hx_index_alloc(k);

  /* ensure the index was allocated. */
  if (!sz)
    throw("failed to allocate %d indices", k);

  /* store the elements of the size index. */
  for (i = 0; i < k; i++)
    sz[i] = x->sz[kx[i]];

  /* allocate the bounding array indices. */
  lower = hx_index_alloc(k);
  upper = hx_index_alloc(k);

  /* ensure the bounding arrays were allocated. */
  if (!lower || !upper)
    throw("failed to allocate bounding arrays");

  /* store the elements of the bounding arrays. */
  for (i = 1; i < k; i++) {
    /* store the upper and lower bound. */
    upper[i] = x->sz[kx[i]] - 1;
    lower[i] = 0;
  }

  /* build a list of packed schedule indices for each slice. */
  ysched = hx_index_scheduled(k - 1, sz + 1, dsched, nsched, sched);

  /* allocate a list of packed schedule indices for the input array. */
  xsched = hx_index_copy(nsched, ysched);

  /* ensure the schedule lists were built successfully. */
  if (!xsched || !ysched)
    throw("failed to build packed schedule arrays");

  /* adjust the elements of the input array's schedule. */
  for (i = 0; i < nsched; i++)
    xsched[i] *= sz[0];

  /* compute the discrete Fourier transform matrix. */
  hx_array_init(&F);
  if (!hx_array_irls_dftmatrix(d, k, nsched, sz, dx, kx, ysched, &F))
    throw("failed to compute dft matrix");

  /* store the dft matrix row and column sizes. */
  n = F.sz[0];
  N = F.sz[1];

  /* allocate a size index for the normal matrix. */
  Asz = hx_index_build(2, n, n);

  /* ensure the index was allocated. */
  if (!Asz)
    throw("failed to allocate matrix size index");

  /* allocate temporary vectors for use during reconstruction. */
  if (!hx_array_alloc(&z, d, 1, &n) ||
      !hx_array_alloc(&y, d, 1, &n) ||
      !hx_array_alloc(&Y, d, 1, &N) ||
      !hx_array_alloc(&w, 0, 1, &N) ||
      !hx_array_alloc(&wx, d, 1, &N) ||
      !hx_array_alloc(&rx, d, 1, &n) ||
      !hx_array_alloc(&A, d, 2, Asz))
    throw("failed to allocate reconstruction arrays");

  /* loop over each reconstruction to be performed. */
  for (is = 0; is < ns; is++) {
    /* store the direct dimension bounds. */
    upper[0] = is;
    lower[0] = is;

    /* slice the indirect dimensions from the input array. */
    if (!hx_array_slice_sched(x, &y, is, nsched, xsched))
      throw("failed to slice sub-matrix %d", is);

    /* reconstruct the current slice. */
    if (!hx_array_irlsfun(&Y, &y, &w, &z, &wx, &rx, &A, niter, pa, pb))
      throw("failed to reconstruct sub-matrix %d", is);

    /* store the reconstructed slice back into the input array. */
    if (!hx_array_store(x, &Y, lower, upper))
      throw("failed to store sub-matrix %d", is);
  }

  /* free the allocated arrays. */
  hx_array_free(&A);
  hx_array_free(&F);
  hx_array_free(&Y);
  hx_array_free(&y);
  hx_array_free(&w);
  hx_array_free(&z);
  hx_array_free(&wx);
  hx_array_free(&rx);

  /* free the allocated indices. */
  hx_index_free(xsched);
  hx_index_free(ysched);
  hx_index_free(lower);
  hx_index_free(upper);
  hx_index_free(Asz);
  hx_index_free(sz);

  /* return success. */
  return 1;
}

