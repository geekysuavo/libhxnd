
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

/* define the number of corrections to store that approximate the inverse
 * of the hessian matrix.
 */
#define HX_MAXENT_NUM_CORRECTIONS  5

/* functional_def: structure definition for looking up entropy functionals
 * by string name or enumerated type.
 */
struct functional_def {
  /* @name: entropy functional string name.
   * @type: entropy functional enumerated type.
   */
  const char *name;
  enum hx_entropy_type type;

  /* @f: entropy functional scalar value function pointer.
   * @df: entropy functional partial derivative function pointer.
   */
  hx_entropy_functional f;
  hx_entropy_functional df;
};

/* functionals: structure array of all supported entropy functionals.
 */
const struct functional_def functionals[] = {
  /* l1-norm. */
  { HX_ENTROPY_NAME_NORM,
    HX_ENTROPY_TYPE_NORM,
    &hx_entropy_norm_f,
    &hx_entropy_norm_df
  },

  /* shannon entropy. */
  { HX_ENTROPY_NAME_SHANNON,
    HX_ENTROPY_TYPE_SHANNON,
    &hx_entropy_shannon_f,
    &hx_entropy_shannon_df
  },

  /* skilling entropy. */
  { HX_ENTROPY_NAME_SKILLING,
    HX_ENTROPY_TYPE_SKILLING,
    &hx_entropy_skilling_f,
    &hx_entropy_skilling_df
  },

  /* hoch/hore entropy. */
  { HX_ENTROPY_NAME_HOCH,
    HX_ENTROPY_TYPE_HOCH,
    &hx_entropy_hoch_f,
    &hx_entropy_hoch_df
  },

  /* null termination. */
  { NULL, HX_ENTROPY_TYPE_UNDEFINED, NULL, NULL }
};

/* hx_entropy_lookup_type(): return the enumerated entropy functional type
 * based on a specified string representation.
 * @name: the entropy functional name string.
 */
enum hx_entropy_type hx_entropy_lookup_type (const char *name) {
  /* declare a required variable. */
  unsigned int i;

  /* return an undefined entropy type if the name is null. */
  if (!name)
    return HX_ENTROPY_TYPE_UNDEFINED;

  /* loop over all supported entropy functionals. */
  for (i = 0; functionals[i].name; i++) {
    /* check if the functional name matches. */
    if (strcmp(name, functionals[i].name) == 0)
      return functionals[i].type;
  }

  /* return an undefined entropy type. */
  return HX_ENTROPY_TYPE_UNDEFINED;
}

/* hx_entropy_get_functionals: get the entropy functional function pointers.
 * @type: entropy functional type to query.
 * @f: location to store the entropy function pointer.
 * @df: location to store the derivative function pointer.
 */
int hx_entropy_get_functionals (enum hx_entropy_type type,
                                hx_entropy_functional *f,
                                hx_entropy_functional *df) {
  /* declare a required variable. */
  unsigned int i;

  /* loop over all supported entropy functional types. */
  for (i = 0; functionals[i].name; i++) {
    /* check if the functional type matches. */
    if (functionals[i].type == type) {
      /* store the function pointers. */
      *f = functionals[i].f;
      *df = functionals[i].df;

      /* return success. */
      return 1;
    }
  }

  /* return failure. */
  return 0;
}

/* hx_entropy_norm_f(): compute the scalar entropy of a hypercomplex scalar
 * using an l1-norm functional.
 */
int hx_entropy_norm_f (real *x, real *S, int n) {
  /* compute the norm of the hypercomplex input. */
  *S = hx_data_real_norm(x, n);

  /* return success. */
  return 1;
}

/* hx_entropy_norm_df(): compute the complex partial derivative of the
 * entropy of a hypercomplex scalar using an l1-norm functional.
 */
int hx_entropy_norm_df (real *x, real *S, int n) {
  /* declare required variables. */
  real Snrm;
  int i;

  /* compute the norm of the hypercomplex input. */
  Snrm = hx_data_real_norm(x, n);

  /* compute the array elements of the complex derivative. */
  for (i = 0; i < n; i++)
    S[i] = x[i] / Snrm;

  /* return success. */
  return 1;
}

/* hx_entropy_shannon_f(): compute the scalar entropy of a hypercomplex scalar
 * using a shannon entropy functional.
 */
int hx_entropy_shannon_f (real *x, real *S, int n) {
  /* declare required variables. */
  real Snrm;

  /* compute the entropy of the hypercomplex input. */
  Snrm = hx_data_real_norm(x, n);
  *S = Snrm * log(Snrm);

  /* return success. */
  return 1;
}

/* hx_entropy_shannon_df(): compute the complex partial derivative of the
 * entropy of a hypercomplex scalar using a shannon entropy functional.
 */
int hx_entropy_shannon_df (real *x, real *S, int n) {
  /* declare required variables. */
  real Snrm;
  int i;

  /* compute the norm of the hypercomplex input. */
  Snrm = hx_data_real_norm(x, n);
  Snrm = (log(Snrm) + 1.0) / Snrm;

  /* compute the array elements of the complex derivative. */
  for (i = 0; i < n; i++)
    S[i] = Snrm * x[i];

  /* return success. */
  return 1;
}

/* hx_entropy_skilling_f(): compute the scalar entropy of a hypercomplex scalar
 * using a skilling entropy functional.
 */
int hx_entropy_skilling_f (real *x, real *S, int n) {
  /* declare required variables. */
  real Snrm;

  /* compute the entropy of the hypercomplex input. */
  Snrm = hx_data_real_norm(x, n);
  *S = Snrm * log(Snrm) - Snrm;

  /* return success. */
  return 1;
}

/* hx_entropy_skilling_df(): compute the complex partial derivative of the
 * entropy of a hypercomplex scalar using a skilling entropy functional.
 */
int hx_entropy_skilling_df (real *x, real *S, int n) {
  /* declare required variables. */
  real Snrm;
  int i;

  /* compute the norm of the hypercomplex input. */
  Snrm = hx_data_real_norm(x, n);
  Snrm = log(Snrm) / Snrm;

  /* compute the array elements of the complex derivative. */
  for (i = 0; i < n; i++)
    S[i] = Snrm * x[i];

  /* return success. */
  return 1;
}

/* hx_entropy_hoch_f(): compute the scalar entropy of a hypercomplex scalar
 * using a hoch/hore spin-half entropy functional.
 */
int hx_entropy_hoch_f (real *x, real *S, int n) {
  /* declare required variables. */
  real Snrm;

  /* compute the norm of the hypercomplex input. */
  Snrm = hx_data_real_norm(x, n);

  /* compute the entropy of the hypercomplex input. */
  *S = Snrm * log(Snrm / 2.0 + sqrt(1.0 + Snrm * Snrm / 4.0))
     - sqrt(4.0 + Snrm * Snrm);

  /* return success. */
  return 1;
}

/* hx_entropy_hoch_df(): compute the complex partial derivative of the
 * entropy of a hypercomplex scalar using a hoch/hore spin-half entropy
 * functional.
 */
int hx_entropy_hoch_df (real *x, real *S, int n) {
  /* declare required variables. */
  real Snrm;
  int i;

  /* compute the norm of the hypercomplex input. */
  Snrm = hx_data_real_norm(x, n);
  Snrm = log(Snrm / 2.0 + sqrt(1.0 + Snrm * Snrm / 4.0)) / Snrm;

  /* compute the array elements of the complex derivative. */
  for (i = 0; i < n; i++)
    S[i] = Snrm * x[i];

  /* return success. */
  return 1;
}

/* hx_array_ffmfn(): compute the fast forward maximum entropy reconstruction
 * of a single submatrix sliced from a hypercomplex array.
 *  - refer to comments within hx_array_ffm() for argument descriptions.
 */
int hx_array_ffmfn (hx_array *X, hx_array *x, hx_array *x0,
                    hx_array *g, hx_array *h, hx_array *S, hx_array *Y,
                    hx_array *alpha, hx_array *beta, hx_array *rho,
                    hx_entropy_functional f,
                    hx_entropy_functional df,
                    int niter) {
  /* FIXME: implement hx_array_ffmfn() */

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
   * @Y: array of dense frequency-domain estimates.
   * @y: array of dense time-domain estimates.
   * @y0: vector of sparse time-domain values.
   * @g: current iteration gradient vector.
   * @h: current iteration step direction.
   * @W: matrix of spectral hessian corrections.
   * @Z: matrix of gradient hessian corrections.
   * @alpha, @beta, @rho: vectors used in hessian update.
   * @sz: size of each slice to be reconstructed.
   * @ysched: array of packed schedule indices.
   * @lower: reconstruction slice lower bound index.
   * @upper: reconstruction slice upper bound index.
   * @f: function pointer for the entropy functional.
   * @df: function pointer for the entropy derivative.
   * @d: number of array/slice algebraic dimensions.
   * @k: number of array/slice topological dimensions.
   * @i: general purpose loop counter.
   * @is: loop counter of reconstruction slices.
   * @ns: number of slices to reconstruct.
   * @M: number of hessian corrections to store.
   * @N: number of frequency-domain data points.
   */
  hx_array Y, y, y0, g, h, W, Z, alpha, beta, rho;
  hx_index Wsz, sz, ysched, lower, upper;
  hx_entropy_functional f, df;
  int d, k, i, is, ns, M, N;

  /* ensure the schedule is allocated and properly sized. */
  if (!sched || dsched < 1 || nsched < 1)
    throw("invalid schedule configuration (%dx%d)", dsched, nsched);

  /* ensure the iteration count is in bounds. */
  if (niter < 1)
    throw("iteration count %d out of bounds [1,inf)", niter);

  /* retrieve the entropy functional function pointers. */
  if (!hx_entropy_get_functionals(type, &f, &df))
    throw("failed to retrieve entropy functionals");

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

  /* ensure the schedule list was built successfully. */
  if (!ysched)
    throw("failed to build packed schedule array");

  /* compute the number of frequency-domain indices. */
  for (i = 1, N = 1; i < k; i++)
    N *= sz[i];

  /* set the number of hessian corrections. */
  M = HX_MAXENT_NUM_CORRECTIONS;

  /* allocate a size index for the hessian corrections. */
  Wsz = hx_index_build(2, N, M);

  /* ensure the index was allocated. */
  if (!Wsz)
    throw("failed to allocate matrix size index");

  /* allocate temporary arrays for use during reconstruction. */
  if (!hx_array_alloc(&Y, d, k, sz) ||
      !hx_array_alloc(&y, d, k, sz) ||
      !hx_array_alloc(&y0, d, 1, &nsched) ||
      !hx_array_alloc(&g, d, 1, &N) ||
      !hx_array_alloc(&h, d, 1, &N) ||
      !hx_array_alloc(&W, d, 2, Wsz) ||
      !hx_array_alloc(&Z, d, 2, Wsz) ||
      !hx_array_alloc(&alpha, d, 1, &M) ||
      !hx_array_alloc(&beta, d, 1, &M) ||
      !hx_array_alloc(&rho, d, 1, &M))
    throw("failed to allocate reconstruction arrays");

  /* loop over each reconstruction to be performed. */
  for (is = 0; is < ns; is++) {
    /* store the direct dimension bounds. */
    upper[0] = is;
    lower[0] = is;

    /* slice the indirect dimensions from the input array. */
    if (!hx_array_slice(x, &y, lower, upper))
      throw("failed to slice sub-matrix %d", is);

    /* build a vector containing only the initially sampled points. */
    if (!hx_array_slice_sched(&y, &y0, 0, nsched, ysched))
      throw("failed to initialize time-domain vector %d", is);

    /* reconstruct the current slice. */
    if (!hx_array_ffmfn(&Y, &y, &y0, &g, &h, &W, &Z,
                        &alpha, &beta, &rho, f, df, niter))
      throw("failed to reconstruct sub-matrix %d", is);

    /* store the reconstructed slice back into the input array. */
    if (!hx_array_store(x, &y, lower, upper))
      throw("failed to store sub-matrix %d", is);
  }

  /* free the allocated arrays. */
  hx_array_free(&Y);
  hx_array_free(&y);
  hx_array_free(&y0);
  hx_array_free(&g);
  hx_array_free(&h);
  hx_array_free(&W);
  hx_array_free(&Z);
  hx_array_free(&alpha);
  hx_array_free(&beta);
  hx_array_free(&rho);

  /* free the allocated indices. */
  hx_index_free(ysched);
  hx_index_free(lower);
  hx_index_free(upper);
  hx_index_free(Wsz);
  hx_index_free(sz);

  /* return success. */
  return 1;
}

