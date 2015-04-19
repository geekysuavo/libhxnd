
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

/* declare a few fixed constants used during optimization.
 */
#define HX_IRLS_EPSILON     0.0001
#define HX_IRLS_LAMBDA_MIN  1.0e-3
#define HX_IRLS_LAMBDA_MAX  1.0e+9

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
    hx_index_unpack(k - 1, sz + 1, idxi + 1, sched[i]);

    /* loop over the output matrix columns. */
    for (j = 0; j < N; j++) {
      /* unpack the current frequency-domain index. */
      hx_index_unpack(k - 1, sz + 1, idxk + 1, j);

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
            !hx_scalar_zero(&ph) ||
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

/* hx_array_irls_dft(): compute the discrete Fourier transform, or its
 * inverse, between a nonuniformly sampled time-domain array and a
 * complete frequency-domain array.
 * @F: matrix of Fourier coefficients.
 * @X: vector of frequency-domain values.
 * @x: vector of time-domain values.
 * @dir: transform direction.
 */
int hx_array_irls_dft (hx_array *F, hx_array *X, hx_array *x, real dir) {
  /* declare a few required variables:
   * @i, @n: transform row count and loop index.
   * @k, @N: transform column count and loop index.
   */
  int i, k, n, N, idxf;

  /* store the sizes of the transform matrix. */
  n = F->sz[0];
  N = F->sz[1];

  /* determine the transform direction. */
  if (dir == HX_FFT_FORWARD) {
    /* zero the output vector elements. */
    hx_array_zero(X);

    /* create a team of threads to compute each inner product. */
    #pragma omp parallel private(i, k, idxf)
    {
      /* declare a few required thread-local variables:
       * @Fh: conjugated copy of the current matrix element.
       */
      hx_scalar Fh;

      /* allocate memory for the conjugated scalar. */
      if (!hx_scalar_alloc(&Fh, F->d))
        raise("failed to allocate temporary %d-scalar", F->d);

      /* loop over the transform inner products. */
      #pragma omp for
      for (k = 0; k < N; k++) {
        /* loop over the inner product terms. */
        for (i = 0; i < n; i++) {
          /* compute the matrix coefficient index. */
          idxf = (i + k * n) * F->n;

          /* copy-semiconjugate the current matrix coefficient. */
          hx_data_semiconj(F->x + idxf, Fh.x, F->d, F->n);

          /* compute the current inner product term. */
          hx_data_mul(Fh.x,
                      x->x + i * F->n,
                      X->x + k * F->n,
                      F->d, F->n, F->tbl);
        }
      }

      /* free the temporary scalar. */
      hx_scalar_free(&Fh);
    }
  }
  else {
    /* zero the output vector elements. */
    hx_array_zero(x);

    /* loop over the transform inner products. */
    #pragma omp parallel for private(i, k, idxf)
    for (i = 0; i < n; i++) {
      /* loop over the inner product terms. */
      for (k = 0; k < N; k++) {
        /* compute the matrix coefficient index. */
        idxf = (i + k * n) * F->n;

        /* compute the current inner product term. */
        hx_data_mul(F->x + idxf,
                    X->x + k * F->n,
                    x->x + i * F->n,
                    F->d, F->n, F->tbl);
      }
    }
  }

  /* return success. */
  return 1;
}

/* hx_array_irls_reweight(): compute the set of weights for a spectral
 * estimate vector.
 * @X: vector of input frequency-domain spectral estimates.
 * @w: vector of output real weighting coefficients.
 * @p: current norm p-value for the weighting.
 */
int hx_array_irls_reweight (hx_array *X, hx_array *w, real p) {
  /* declare a few required variables:
   * @Xnrm: norm of the current spectral point.
   * @i: weight vector index, loop counter.
   */
  real Xnrm;
  int i;

  /* loop over the vector elements. */
  for (i = 0; i < w->len; i++) {
    /* compute the norm of the current spectral point. */
    Xnrm = hx_data_real_norm(X->x + i * X->n, X->n);

    /* compute the new weighting value. */
    w->x[i] = 1.0 / (pow(Xnrm, 2.0 - p) + HX_IRLS_EPSILON);
  }

  /* return success. */
  return 1;
}

/* hx_array_irls_sumsq(): compute the sums of squares of the weighted
 * spectral estimate and the residuals, and adjust the weights accordingly.
 * @X: vector of frequency-domain estimates.
 * @x: vector of time-domain acquired values.
 * @z: vector of time-domain estimated values.
 * @w: vector of real weights.
 */
int hx_array_irls_sumsq (hx_array *X, hx_array *x,
                         hx_array *z, hx_array *w) {
  /* declare a few required variables:
   * @rx: sum of squares of the residuals.
   * @wx: sum of squares of the weighted spectrum.
   * @i: vector index, loop counter.
   */
  real rx, wx, rxi, wxi, lambda;
  hx_scalar r;
  int i;

  /* allocate a hypercomplex scalar to hold residuals. */
  if (!hx_scalar_alloc(&r, x->d))
    throw("failed to allocate %d-scalar", x->d);

  /* loop over the first sum of squares computation. */
  for (i = 0, rx = 0.0; i < x->len; i += x->n) {
    /* compute the residual term. */
    hx_data_add(x->x + i, z->x + i, r.x, -1.0, x->d, x->n);

    /* update the sum of squares. */
    rxi = hx_data_real_norm(r.x, r.n);
    rx += rxi * rxi;
  }

  /* loop over the second sum of squares computation. */
  for (i = 0, wx = 0.0; i < w->len; i++) {
    /* update the sum of squares. */
    wxi = w->x[i] * hx_data_real_norm(X->x + i * X->n, X->n);
    wx += wxi * wxi;
  }

  /* compute the lagrange multiplier. */
  lambda = rx / wx;

  /* ensure lambda is never too small. */
  if (lambda < HX_IRLS_LAMBDA_MIN)
    lambda = HX_IRLS_LAMBDA_MIN;

  /* ensure lambda is never too large. */
  if (lambda > HX_IRLS_LAMBDA_MAX)
    lambda = HX_IRLS_LAMBDA_MAX;

  /* scale the vector of weights by the lagrange multiplier and then
   * invert the vector element such that the weight vector contains
   * diag(inv(W)).
   */
  for (i = 0; i < w->len; i++)
    w->x[i] = 1.0 / (lambda * w->x[i]);

  /* free the temporary scalar. */
  hx_scalar_free(&r);

  /* return success. */
  return 1;
}

/* hx_array_irls_gramian(): compute the regression coefficient matrix which
 * will be decomposed during solution for the new time-domain estimate.
 * @F: matrix of Fourier coefficients.
 * @w: vector of weighting factors.
 * @A: output gram matrix.
 */
int hx_array_irls_gramian (hx_array *F, hx_array *w, hx_array *A) {
  /* declare a few required variables:
   * @n: number of transform matrix rows.
   * @N: number of transform matrix columns.
   */
  int n, N;

  /* store the sizes of the transform matrix. */
  n = F->sz[0];
  N = F->sz[1];

  /* zero the output matrix elements. */
  hx_array_zero(A);

  /* create a team of threads to compute matrix elements in parallel. */
  #pragma omp parallel
  {
    /* declare a few required thread-local variables:
     * @i, @j, @k: matrix-matrix product loop counters.
     * @idx: coefficient index of the unconjugated transform scalar.
     * @idxh: coefficient index of the conjugated transform scalar.
     * @idxa: coefficient index of the output matrix element.
     * @Fh: temporary scalar holding the conjugated, weighted element.
     */
    int i, j, k, idx, idxh, idxa;
    hx_scalar Fh;

    /* allocate a hypercomplex scalar to hold residuals. */
    if (!hx_scalar_alloc(&Fh, F->d))
      raise("failed to allocate %d-scalar", F->d);

    /* loop over the design matrix columns. */
    #pragma omp for
    for (j = 0; j < n; j++) {
      /* loop over the set of inner products. */
      for (k = 0; k < N; k++) {
        /* compute the conjugated matrix element index. */
        idxh = (j + k * n) * F->n;

        /* copy-semiconjugate the current transform matrix element. */
        hx_data_semiconj(F->x + idxh, Fh.x, F->d, F->n);

        /* scale the conjugated transform matrix element. */
        hx_data_add(NULL, Fh.x, Fh.x, w->x[k], Fh.d, Fh.n);

        /* loop over the design matrix rows. */
        for (i = 0; i < n; i++) {
          /* compute the output matrix element index. */
          idxa = (i + j * n) * A->n;

          /* compute the unconjugated matrix element index. */
          idx = (i + k * n) * F->n;

          /* compute the current inner product term. */
          hx_data_mul(F->x + idx, Fh.x, A->x + idxa,
                      F->d, F->n, F->tbl);
        }
      }

      /* tack on the identity matrix. */
      A->x[(j + j * n) * A->n] += 1.0;
    }

    /* free the temporary scalar. */
    hx_scalar_free(&Fh);
  }

  /* return success. */
  return 1;
}

/* hx_array_irls_solve(): solve the linear system of equations,
 *   A * x = L * L^T * x = b
 */
int hx_array_irls_solve (hx_array *A, hx_array *x, hx_array *b) {
  /* declare a few required variables:
   */
  int i, j, k, n, idxjj;
  hx_scalar tmp, Lh;
  real Anrm, diag;

  /* store the problem size locally. */
  n = A->sz[0];

  /* allocate the temporary substitution scalar. */
  if (!hx_scalar_alloc(&tmp, A->d) ||
      !hx_scalar_alloc(&Lh, A->d))
    throw("failed to allocate temporary %d-scalars", A->d);

  /* loop serially over the rows of the matrix. */
  for (j = idxjj = 0; j < n; j++, idxjj += A->n * (n + 1)) {
    /* initialize the value of the new diagonal. */
    diag = A->x[idxjj];

    /* loop to compute the new diagonal element. */
    for (k = 0; k < j; k++) {
      /* compute the norm of the preceeding off-diagonals. */
      Anrm = hx_data_real_norm(A->x + (j + k * n) * A->n, A->n);
      diag -= Anrm * Anrm;
    }

    /* ensure the new diagonal element is positive. */
    if (diag <= 0.0)
      throw("pivot %d is %.3le", j, diag);

    /* store the newly computed diagonal element. */
    hx_data_zero(A->x + idxjj, A->n);
    A->x[idxjj] = sqrt(diag);
    diag = A->x[idxjj];

    /* create a team of threads to execute the inner loop. */
    #pragma omp parallel private(i, k)
    {
      /* declare a few required variables:
       * @sum: temporary hypercomplex summand.
       * @Ah: conjugated matrix element.
       */
      int idxij, idxik, idxjk;
      hx_scalar sum, Ah;

      /* allocate the temporary scalars. */
      if (!hx_scalar_alloc(&sum, A->d) ||
          !hx_scalar_alloc(&Ah, A->d))
        raise("failed to allocate %d-scalars", A->d);

      /* loop over the elements of the column update. */
      #pragma omp for
      for (i = j + 1; i < n; i++) {
        /* build the output matrix element index. */
        idxij = (i + j * n) * A->n;

        /* initialize the sum of the off-diagonals. */
        hx_scalar_zero(&sum);

        /* compute the sum of the off-diagonals. */
        for (k = 0; k < j; k++) {
          /* build the two required matrix coefficient indices. */
          idxik = (i + k * n) * A->n;
          idxjk = (j + k * n) * A->n;

          /* copy-conjugate the second matrix element. */
          hx_data_semiconj(A->x + idxjk, Ah.x, Ah.d, Ah.n);

          /* compute the current product term. */
          hx_data_mul(A->x + idxik, Ah.x, sum.x, sum.d, sum.n, sum.tbl);
        }

        /* compute the final matrix element. */
        hx_data_add(A->x + idxij, sum.x, A->x + idxij, -1.0, A->d, A->n);
        hx_data_add(NULL, A->x + idxij, A->x + idxij, 1.0 / diag, A->d, A->n);
      }

      /* free the allocated scalars. */
      hx_scalar_free(&sum);
      hx_scalar_free(&Ah);
    }
  }

  /* initialize the output vector. */
  hx_array_zero(x);

  /* perform the first step of substitution. */
  for (j = idxjj = 0; j < n; j++, idxjj += A->n * (n + 1)) {
    /* initialize the temporary sum. */
    hx_scalar_zero(&tmp);

    /* compute the temporary sum. */
    for (k = j - 1; k >= 0; k--) {
      /* compute the current term. */
      hx_data_mul(A->x + (j + k * n) * A->n, x->x + k * x->n, tmp.x,
                  A->d, A->n, A->tbl);
    }

    /* store the final result. */
    hx_data_add(b->x + j * b->n, tmp.x, x->x + j * x->n, -1.0, A->d, A->n);
    hx_data_add(NULL, x->x + j * x->n, x->x + j * x->n,
                1.0 / A->x[idxjj], A->d, A->n);
  }

  /* perform the second step of substitution. */
  for (j = n - 1, idxjj = (n * n - 1) * A->n; j >= 0;
       j--, idxjj -= A->n * (n + 1)) {
    /* initialize the temporary sum. */
    hx_scalar_zero(&tmp);

    /* compute the temporary sum. */
    for (k = j + 1; k < n; k++) {
      /* copy-semiconjugate the current matrix element. */
      hx_data_semiconj(A->x + (k + j * n) * A->n, Lh.x, A->d, A->n);

      /* compute the current term. */
      hx_data_mul(Lh.x, x->x + k * x->n, tmp.x,
                  A->d, A->n, A->tbl);
    }

    /* store the final result. */
    hx_data_add(x->x + j * x->n, tmp.x, x->x + j * x->n, -1.0, A->d, A->n);
    hx_data_add(NULL, x->x + j * x->n, x->x + j * x->n,
                1.0 / A->x[idxjj], A->d, A->n);
  }

  /* free the temporary scalars. */
  hx_scalar_free(&tmp);
  hx_scalar_free(&Lh);

  /* return success. */
  return 1;
}

/* hx_array_irls_weight_spect(): weight the spectral estimate produced by
 * fourier transformation of the regression solution.
 * @X: current frequency-domain spectral estimate.
 * @w: current set of (inverted) weights.
 */
int hx_array_irls_weight_spect (hx_array *X, hx_array *w) {
  /* declare a required variable. */
  int i, idx;

  /* loop over the vector elements. */
  for (i = idx = 0; i < w->len; i++, idx += X->n)
    hx_data_add(NULL, X->x + idx, X->x + idx, w->x[i], X->d, X->n);

  /* return success. */
  return 1;
}

/* hx_array_irlsfn(): compute a single irls reconstruction.
 * @F: matrix of Fourier coefficients.
 * @X: vector of frequency-domain spectral estimates.
 * @x: vector of time-domain acquired values.
 * @w: vector of real spectral weights.
 * @z: vector of time-domain estimates.
 * @A: hermitian positive definite array of regression factors.
 * @niter, @pa, @pb: see hx_array_irls().
 */
int hx_array_irlsfn (hx_array *F, hx_array *X, hx_array *x,
                     hx_array *w, hx_array *z, hx_array *A,
                     int niter, real pa, real pb) {
  /* declare a few required variables:
   * @p: current iteration norm p-value.
   * @dp: change in p-value per iteration.
   * @iiter: iteration loop counter.
   */
  real p, dp;
  int iiter;

  /* compute the change in norm p-value per iteration. */
  dp = (pb - pa) / (real) niter;

  /* compute the initial spectral estimate. */
  if (!hx_array_irls_dft(F, X, x, HX_FFT_FORWARD))
    throw("failed to compute initial dft");

  /* store the initial time-domain estimate. */
  if (!hx_data_copy(x->x, z->x, x->len))
    throw("failed to initialize time-domain vector");

  /* loop over the reconstruction iterations. */
  for (iiter = 0; iiter < niter; iiter++) {
    /* compute the current norm p-value. */
    p = pa + (real) iiter * dp;

    /* compute the current weights. */
    if (!hx_array_irls_reweight(X, w, p))
      throw("failed to compute new weights");

    /* adjust the weights to equalize the terms of the minimization. */
    if (!hx_array_irls_sumsq(X, x, z, w))
      throw("failed to adjust new weights");

    /* compute the regression coefficient matrix. */
    if (!hx_array_irls_gramian(F, w, A))
      throw("failed to compute gram matrix");

    /* decompose the design matrix and solve for the time-domain vector. */
    if (!hx_array_irls_solve(A, z, x))
      throw("failed to solve linear system");

    /* compute the unweighted frequency-domain estimate. */
    if (!hx_array_irls_dft(F, X, z, HX_FFT_FORWARD))
      throw("failed to compute forward dft");

    /* weight the spectral estimate. */
    if (!hx_array_irls_weight_spect(X, w))
      throw("failed to weight spectral estimate");

    /* compute the new time-domain estimate. */
    if (!hx_array_irls_dft(F, X, z, HX_FFT_REVERSE))
      throw("failed to compute inverse dft");
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
   * @Asz: size index of the regression coefficient matrix.
   * @sz: size of each slice to be reconstructed.
   * @xsched: array of packed schedule indices within @x.
   * @ysched: array of packed schedule indices within @y.
   * @lower: reconstruction slice lower bound index.
   * @upper: reconstruction slice upper bound index.
   * @A: regression coefficient matrix.
   * @F: discrete Fourier transform matrix.
   * @Y: spectral estimate vector.
   * @y: time-domain data vector.
   * @w: real spectral weight vector.
   * @z: time-domain estimate vector.
   * @d: number of array/slice algebraic dimensions.
   * @k: number of array/slice topological dimensions.
   * @i: general purpose loop counter.
   * @is: loop counter of reconstruction slices.
   * @ns: number of slices to reconstruct.
   * @n: total number of sampled time-domain points per slice.
   * @N: total number of frequency-domain points per slice.
   */
  hx_index Asz, sz, xsched, ysched, lower, upper;
  hx_array A, F, Y, y, w, z;
  int d, k, i, is, ns, n, N;

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
    if (!hx_array_irlsfn(&F, &Y, &y, &w, &z, &A, niter, pa, pb))
      throw("failed to reconstruct sub-matrix %d", is);

    /* store the reconstructed slice back into the input array. */
    if (!hx_array_store(x, &Y, lower, upper))
      throw("failed to store sub-matrix %d", is);
  }

  /* inverse fourier transform the shifted result. */
  for (i = 1; i < k; i++) {
    /* fourier transform the current dimension. */
    if (!hx_array_ifft(x, dx[i], kx[i]))
      throw("failed to apply final inverse fft");
  }

  /* free the allocated arrays. */
  hx_array_free(&A);
  hx_array_free(&F);
  hx_array_free(&Y);
  hx_array_free(&y);
  hx_array_free(&w);
  hx_array_free(&z);

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

