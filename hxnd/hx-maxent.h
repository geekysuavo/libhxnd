
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

/* ensure once-only inclusion. */
#ifndef __HXND_HX_MAXENT_H__
#define __HXND_HX_MAXENT_H__

/* define string constants for supported entropy functionals.
 */
#define HX_ENTROPY_NAME_NORM      "norm"
#define HX_ENTROPY_NAME_SHANNON   "shannon"
#define HX_ENTROPY_NAME_SKILLING  "skilling"
#define HX_ENTROPY_NAME_HOCH      "hoch"

/* hx_entropy_type: enumerated type for entropy functionals.
 */
enum hx_entropy_type {
  HX_ENTROPY_TYPE_UNDEFINED,
  HX_ENTROPY_TYPE_NORM,
  HX_ENTROPY_TYPE_SHANNON,
  HX_ENTROPY_TYPE_SKILLING,
  HX_ENTROPY_TYPE_HOCH
};

/* hx_entropy_functional: function prototype for all entropy functionals.
 * @x: hypercomplex spectral intensity array raw data.
 * @S: output entropy or entropy derivative pointer.
 * @n: number of coefficients per hypercomplex scalar.
 */
typedef int (*hx_entropy_functional) (real *x, real *S, int n);

/* function declarations: */

enum hx_entropy_type hx_entropy_lookup_type (const char *name);

int hx_entropy_norm_f (real *x, real *S, int n);

int hx_entropy_norm_df (real *x, real *S, int n);

int hx_entropy_shannon_f (real *x, real *S, int n);

int hx_entropy_shannon_df (real *x, real *S, int n);

int hx_entropy_skilling_f (real *x, real *S, int n);

int hx_entropy_skilling_df (real *x, real *S, int n);

int hx_entropy_hoch_f (real *x, real *S, int n);

int hx_entropy_hoch_df (real *x, real *S, int n);

int hx_array_ffm (hx_array *x, hx_index dx, hx_index kx,
                  int dsched, int nsched, hx_index sched,
                  int niter, enum hx_entropy_type type);

#endif /* __HXND_HX_MAXENT_H__ */

