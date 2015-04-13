
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
  for (i = 0, N = 1; i < k; i++)
    N *= sz[i];

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
    hx_index_unpack(k, sz, idxi, sched[i]);

    /* loop over the output matrix columns. */
    for (j = 0; j < N; j++) {
      /* unpack the current frequency-domain index. */
      hx_index_unpack(k, sz, idxk, j);

      /* initialize the final scalar phase factor. */
      hx_scalar_zero(&ph);
      ph.x[0] = 1.0;

      /* loop over the transform dimensions. */
      for (di = 0; di < d; di++) {
        /* compute the current dimension phase angle. */
        theta = 2.0 * M_PI;
        theta *= (real) idxi[di];
        theta *= (real) idxk[di];
        theta /= (real) sz[di];

        /* 1. compute the current dimension phase factor.
         * 2. copy the current phase factor into a temporary scalar.
         * 3. multiply the phase factor into the final phase factor.
         */
        if (!hx_data_copy(ph.x, phtmp.x, ph.n) ||
            !hx_scalar_phasor(&phd, di, theta) ||
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

  /* FIXME: implement hx_array_irls() */

  /* return success. */
  return 1;
}

