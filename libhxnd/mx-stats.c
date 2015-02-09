
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
  for (i = 1, val[0] = x->x[0]; i < x->len; i++) {
    /* check if the current value is a new minimum. */
    if (x->x[i] < val[0])
      val[0] = x->x[i];
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
  for (i = 1, val[0] = x->x[0]; i < x->len; i++) {
    /* check if the current value is a new maximum. */
    if (x->x[i] > val[0])
      val[0] = x->x[i];
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
  val[0] = max - min;
  return 1;
}

/* mx_stats_med_qsel(): quickselect algorithm for locating the
 * @k-th order statistic from an array.
 */
real mx_stats_med_qsel (hx_array *x, int k) {
  /* declare a few required variables:
   * @il, @ir: left and right array indices.
   * @ipiv: pivot array index.
   */
  int i, il, ir, ipiv;
  real xpiv, swp;

  /* initialize the array indices. */
  il = 0;
  ir = x->len - 1;

  /* loop until the left and right indices converge. */
  while (il != ir) {
    /* compute the pivot index. this is less ideal than using a
     * pseudorandom pivot, but it's better than nothing.
     */
    ipiv = il + (ir - il) / 2;
    xpiv = x->x[ipiv];

    /* swap the pivot value to the end. */
    x->x[ipiv] = x->x[ir];
    x->x[ir] = xpiv;

    /* partition the sub-array. */
    for (i = ipiv = il; i < ir; i++) {
      /* swap all values less than the pivot. */
      if (x->x[i] < xpiv) {
        /* swap the store index and the current one. */
        swp = x->x[ipiv];
        x->x[ipiv] = x->x[i];
        x->x[i] = swp;

        /* increment the store index. */
        ipiv++;
      }
    }

    /* swap the pivot to its final location. */
    swp = x->x[ir];
    x->x[ir] = x->x[ipiv];
    x->x[ipiv] = swp;

    /* recompute the left and right indices. */
    if (ipiv == k)
      return x->x[k];
    else if (ipiv > k)
      ir = ipiv - 1;
    else
      il = ipiv + 1;
  }

  /* return the final value. */
  return x->x[il];
}

/* mx_stats_med_fn(): projector callback for mx_stats_med().
 * see hx_array_projector_cb() for more details.
 */
int mx_stats_med_fn (hx_array *x, real *val) {
  /* compute the median using quickselect. */
  if (x->len % 2) {
    /* odd coefficient count. */
    val[0] = mx_stats_med_qsel(x, x->len / 2);
  }
  else {
    /* even coefficient count. */
    val[0] = mx_stats_med_qsel(x, x->len / 2);
    val[0] += mx_stats_med_qsel(x, x->len / 2 - 1);
    val[0] /= 2.0;
  }

  /* return success. */
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
  val[0] = sum / (real) x->len;
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
  val[0] = s / (real) (x->len - 1);
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
  val[0] = sqrt(val[0]);
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
 * see mx_stats_min() for more details.
 */
int mx_stats_max (hx_array *x, int k, hx_array *m) {
  /* execute the maximum projector function. */
  hx_array_init(m);
  return hx_array_projector(x, k, &mx_stats_max_fn, m);
}

/* mx_stats_range(): compute the ranges along a given dimension of
 * a real multidimensional array.
 * see mx_stats_min() for more details.
 */
int mx_stats_range (hx_array *x, int k, hx_array *m) {
  /* execute the maximum projector function. */
  hx_array_init(m);
  return hx_array_projector(x, k, &mx_stats_range_fn, m);
}

/* mx_stats_med(): compute the median values along a given dimension of
 * a real multidimensional array.
 * see mx_stats_min() for more details.
 */
int mx_stats_med (hx_array *x, int k, hx_array *m) {
  /* execute the mean projector function. */
  hx_array_init(m);
  return hx_array_projector(x, k, &mx_stats_med_fn, m);
}

/* mx_stats_mean(): compute the mean values along a given dimension of
 * a real multidimensional array.
 * see mx_stats_min() for more details.
 */
int mx_stats_mean (hx_array *x, int k, hx_array *m) {
  /* execute the mean projector function. */
  hx_array_init(m);
  return hx_array_projector(x, k, &mx_stats_mean_fn, m);
}

/* mx_stats_var(): compute the sample variances along a given dimension of
 * a real multidimensional array.
 * see mx_stats_min() for more details.
 */
int mx_stats_var (hx_array *x, int k, hx_array *m) {
  /* execute the mean projector function. */
  hx_array_init(m);
  return hx_array_projector(x, k, &mx_stats_var_fn, m);
}

/* mx_stats_stdev(): compute the sample standard deviations along a
 * given dimension of a real multidimensional array.
 * see mx_stats_min() for more details.
 */
int mx_stats_stdev (hx_array *x, int k, hx_array *m) {
  /* execute the mean projector function. */
  hx_array_init(m);
  return hx_array_projector(x, k, &mx_stats_stdev_fn, m);
}

