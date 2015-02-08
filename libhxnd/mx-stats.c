
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

/* include the univariate statistics header. */
#include <hxnd/mx-stats.h>

/* mx_stats_min_fn(): projector callback for mx_stats_min().
 * see hx_array_projector_cb() for more details.
 */
int mx_stats_min_fn (hx_array *x, real *val) {
  /* declare a few required variables. */
  int i;

  /* loop over the array elements. */
  for (i = 1, *val = x->x[0]; i < x->len; i++) {
    /* check if the current value is a new minimum. */
    if (x->x[i] < *val)
      *val = x->x[i];
  }

  /* return success. */
  return 1;
}

/* mx_stats_max_fn(): projector callback for mx_stats_max().
 * see hx_array_projector_cb() for more details.
 */
int mx_stats_max_fn (hx_array *x, real *val) {
  /* declare a few required variables. */
  int i;

  /* loop over the array elements. */
  for (i = 1, *val = x->x[0]; i < x->len; i++) {
    /* check if the current value is a new maximum. */
    if (x->x[i] > *val)
      *val = x->x[i];
  }

  /* return success. */
  return 1;
}

/* mx_stats_range_fn(): projector callback for mx_stats_range().
 * see hx_array_projector_cb() for more details.
 */
int mx_stats_range_fn (hx_array *x, real *val) {
  /* declare a few required variables. */
  real min, max;
  int i;

  /* loop over the array elements. */
  for (i = 1, min = max = x->x[0]; i < x->len; i++) {
    /* check if the current value is a new minimum. */
    if (x->x[i] < min)
      min = x->x[i];

    /* check if the current value is a new maximum. */
    if (x->x[i] > max)
      max = x->x[i];
  }

  /* store the result and return success. */
  *val = max - min;
  return 1;
}

/* mx_stats_mean_fn(): projector callback for mx_stats_mean().
 * see hx_array_projector_cb() for more details.
 */
int mx_stats_mean_fn (hx_array *x, real *val) {
  /* declare a few required variables. */
  real sum;
  int i;

  /* sum over the array elements. */
  for (i = 0, sum = 0.0; i < x->len; i++)
    sum += x->x[i];

  /* store the result and return success. */
  *val = sum / (real) x->len;
  return 1;
}

/* mx_stats_var_fn(): projector callback for mx_stats_var().
 * see hx_array_projector_cb() for more details.
 */
int mx_stats_var_fn (hx_array *x, real *val) {
  /* declare a few required variables. */
  real xi, si, m, s;
  int i;

  /* loop over the array elements. */
  for (i = 1, m = x->x[0], s = 0.0; i < x->len; i++) {
    /* extract the currently indexed array element. */
    xi = x->x[i];

    /* compute the new running mean and sum of squares. */
    si = xi - m;
    m += si / (real) (i + 1);
    s += si * (xi - m);
  }

  /* store the result and return success. */
  *val = s / (real) (x->len - 1);
  return 1;
}

/* mx_stats_stdev_fn(): projector callback for mx_stats_stdev().
 * see hx_array_projector_cb() for more details.
 */
int mx_stats_stdev_fn (hx_array *x, real *val) {
  /* compute the variance. */
  if (!mx_stats_var_fn(x, val))
    return 0;

  /* compute the final result and return success. */
  *val = sqrt(*val);
  return 1;
}

/* mx_stats_min(): compute the minimum values along a given dimension of
 * a real multidimensional array.
 * @x: the input array to compute statistics from.
 * @k: the dimension along which to compute the statistic.
 * @m: pointer to the output array.
 */
int mx_stats_min (hx_array *x, int k, hx_array *m) {
  /* execute the minimum projector function. */
  hx_array_init(m);
  return hx_array_projector(x, k, &mx_stats_min_fn, m);
}

/* mx_stats_max(): compute the maximum values along a given dimension of
 * a real multidimensional array.
 * @x: the input array to compute statistics from.
 * @k: the dimension along which to compute the statistic.
 * @m: pointer to the output array.
 */
int mx_stats_max (hx_array *x, int k, hx_array *m) {
  /* execute the maximum projector function. */
  hx_array_init(m);
  return hx_array_projector(x, k, &mx_stats_max_fn, m);
}

/* mx_stats_range(): compute the ranges along a given dimension of
 * a real multidimensional array.
 * @x: the input array to compute statistics from.
 * @k: the dimension along which to compute the statistic.
 * @m: pointer to the output array.
 */
int mx_stats_range (hx_array *x, int k, hx_array *m) {
  /* execute the maximum projector function. */
  hx_array_init(m);
  return hx_array_projector(x, k, &mx_stats_range_fn, m);
}

/* mx_stats_mean(): compute the mean values along a given dimension of
 * a real multidimensional array.
 * @x: the input array to compute statistics from.
 * @k: the dimension along which to compute the statistic.
 * @m: pointer to the output array.
 */
int mx_stats_mean (hx_array *x, int k, hx_array *m) {
  /* execute the mean projector function. */
  hx_array_init(m);
  return hx_array_projector(x, k, &mx_stats_mean_fn, m);
}

/* mx_stats_var(): compute the sample variances along a given dimension of
 * a real multidimensional array.
 * @x: the input array to compute statistics from.
 * @k: the dimension along which to compute the statistic.
 * @m: pointer to the output array.
 */
int mx_stats_var (hx_array *x, int k, hx_array *m) {
  /* execute the mean projector function. */
  hx_array_init(m);
  return hx_array_projector(x, k, &mx_stats_var_fn, m);
}

/* mx_stats_stdev(): compute the sample standard deviations along a
 * given dimension of a real multidimensional array.
 * @x: the input array to compute statistics from.
 * @k: the dimension along which to compute the statistic.
 * @m: pointer to the output array.
 */
int mx_stats_stdev (hx_array *x, int k, hx_array *m) {
  /* execute the mean projector function. */
  hx_array_init(m);
  return hx_array_projector(x, k, &mx_stats_stdev_fn, m);
}

