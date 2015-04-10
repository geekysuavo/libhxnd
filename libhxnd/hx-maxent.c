
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

