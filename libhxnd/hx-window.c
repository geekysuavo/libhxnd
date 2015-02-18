
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

/* include the window function header. */
#include <hxnd/hx-window.h>

/* wnd_def: structure definition for looking up window function types by
 * string name.
 */
struct wnd_def {
  /* @name: window function string name.
   * @type: window function type.
   */
  const char *name;
  enum hx_window_type type;
};

/* windows: structure array of all supported window functions.
 */
const struct wnd_def windows[] = {
  { HX_WINDOW_NAME_SINE,   HX_WINDOW_TYPE_SINE },
  { HX_WINDOW_NAME_EXP,    HX_WINDOW_TYPE_EXP },
  { HX_WINDOW_NAME_GAUSS,  HX_WINDOW_TYPE_GAUSS },
  { HX_WINDOW_NAME_TRAP,   HX_WINDOW_TYPE_TRAP },
  { HX_WINDOW_NAME_TRI,    HX_WINDOW_TYPE_TRI },
  { HX_WINDOW_NAME_BLACK,  HX_WINDOW_TYPE_BLACK },
  { NULL,                  HX_WINDOW_TYPE_UNDEFINED }
};

/* hx_window_lookup_type(): returns the enumerated type based on its string
 * representation.
 * @name: the window name string.
 */
enum hx_window_type hx_window_lookup_type (const char *name) {
  /* declare a required variable. */
  unsigned int i;

  /* loop over all supported window functions. */
  for (i = 0; windows[i].name; i++) {
    /* check if the window function name matches. */
    if (strcmp(name, windows[i].name) == 0)
      return windows[i].type;
  }

  /* return an undefined window type. */
  return HX_WINDOW_TYPE_UNDEFINED;
}

/* hx_window_alloc(): allocates memory for a hypercomplex window function
 * and stores the computed values in the allocated array.
 * @wnd: hypercomplex array pointer for the result.
 * @d: dimensionality for the array.
 * @len: length of the window vector.
 */
int hx_window_alloc (hx_array *wnd, int d, int len) {
  /* store a local copy of the length argument. */
  int sz = len;

  /* allocate memory for the destination array, but only if the output
   * array does not match the configuration needed.
   */
  if ((wnd->d != d || wnd->k != 1 || wnd->sz[0] != len) &&
      !hx_array_alloc(wnd, d, 1, &sz))
    throw("failed to allocate destination array");

  /* return success. */
  return 1;
}

/* hx_window_sine(): allocate and construct a sinusoidal window function.
 * @start: sine starting phase.
 * @end: sine ending phase.
 * @order: sine exponent.
 */
int hx_window_sine (hx_array *wnd, int d, int len, real width,
                    real start, real end, real order) {
  /* declare required variables:
   * @i: integer point index.
   * @fi: fractional point index.
   */
  real fi, xi;
  int i;

  /* check that the start value is within bounds. */
  if (start < 0.0 || start > 1.0)
    throw("start argument %.3f out of bounds [0,1]", start);

  /* check that the end value is within bounds. */
  if (end < 0.0 || end > 1.0)
    throw("end argument %.3f out of bounds [0,1]", end);

  /* check that the order value is within bounds. */
  if (order < 1.0)
    throw("order argument %.3f out of bounds [1,inf)", order);

  /* ensure that the destination array is allocated. */
  if (!hx_window_alloc(wnd, d, len))
    throw("failed to allocate sinusoidal array");

  /* compute the values of the window. */
  for (i = 0; i < len; i++) {
    /* compute the fractional index. */
    fi = ((real) i) / ((real) (len - 1));

    /* compute the current window value. */
    xi = M_PI * (start + (end - start) * fi);
    xi = pow(sin(xi), order);

    /* store the computed window value. */
    wnd->x[i * wnd->n] = xi;
  }

  /* return success. */
  return 1;
}

/* hx_window_exp(): allocate and construct an exponential window function.
 * @lb: exponential line-broadening in Hertz.
 */
int hx_window_exp (hx_array *wnd, int d, int len, real width,
                   real lb) {
  /* declare required variables:
   * @i: integer point index.
   * @t: time index.
   */
  real t, xi;
  int i;

  /* ensure that the destination array is allocated. */
  if (!hx_window_alloc(wnd, d, len))
    throw("failed to allocate exponential array");

  /* compute the values of the window. */
  for (i = 0; i < len; i++) {
    /* compute the time index. */
    t = ((real) i) / width;

    /* compute the current window value. */
    xi = -M_PI * t * lb;
    xi = exp(xi);

    /* store the computed window value. */
    wnd->x[i * wnd->n] = xi;
  }

  /* return success. */
  return 1;
}

/* hx_window_gauss(): allocate and construct a gaussian window function.
 * @invlb: exponential inverse line-broadening in Hertz.
 * @lb: gaussian line-broadening in Hertz.
 * @center: position of gauss maximum.
 */
int hx_window_gauss (hx_array *wnd, int d, int len, real width,
                     real invlb, real lb, real center) {
  /* declare required variables:
   * @i: integer point index.
   * @t0: center time index.
   * @t: time index.
   */
  real t, t0, xi;
  int i;

  /* check that the center value is within bounds. */
  if (center < 0.0 || center > 1.0)
    throw("center argument %.3f out of bounds [0,1]", center);

  /* ensure that the destination array is allocated. */
  if (!hx_window_alloc(wnd, d, len))
    throw("failed to allocate gaussian array");

  /* compute the scaled center value. */
  t0 = center * ((real) (len - 1)) / width;

  /* compute the values of the window. */
  for (i = 0; i < len; i++) {
    /* compute the time index. */
    t = ((real) i) / width;

    /* compute the current window value. */
    xi = M_PI * t * invlb - pow(0.6 * M_PI * lb * (t0 - t), 2.0);
    xi = exp(xi);

    /* store the computed window value. */
    wnd->x[i * wnd->n] = xi;
  }

  /* return success. */
  return 1;
}

/* hx_window_trap(): allocate and construct a trapezoidal window function.
 * @start: size of starting ramp.
 * @end: size of ending ramp.
 */
int hx_window_trap (hx_array *wnd, int d, int len, real width,
                    real start, real end) {
  /* declare required variables:
   * @i: integer point index.
   * @fi: fractional point index.
   */
  real fi, xi;
  int i;

  /* check that the start value is within bounds. */
  if (start < 0.0 || start > 1.0)
    throw("start argument %.3f out of bounds [0,1]", start);

  /* check that the end value is within bounds. */
  if (end < 0.0 || end > 1.0)
    throw("end argument %.3f out of bounds [0,1]", end);

  /* check that the start and end do not cross. */
  if (start > end)
    throw("start argument may not exceed end argument");

  /* ensure that the destination array is allocated. */
  if (!hx_window_alloc(wnd, d, len))
    throw("failed to allocate trapezoidal array");

  /* compute the values of the window. */
  for (i = 0; i < len; i++) {
    /* compute the fractional index. */
    fi = ((real) i) / ((real) (len - 1));

    /* compute the current window value. */
    if (fi >= end) {
      /* final: download slope. */
      xi = fi / (end - 1.0);
    }
    else if (fi >= start) {
      /* middle: plateau. */
      xi = 1.0;
    }
    else {
      /* initial: upward slope. */
      xi = fi / start;
    }

    /* store the computed window value. */
    wnd->x[i * wnd->n] = xi;
  }

  /* return success. */
  return 1;
}

/* hx_window_tri(): allocate and construct a triangular window function.
 * @locus: position of triangle maximum.
 * @start: height of starting point.
 * @end: height of ending point.
 */
int hx_window_tri (hx_array *wnd, int d, int len, real width,
                   real center, real start, real end) {
  /* declare required variables:
   * @i: integer point index.
   * @fi: fractional point index.
   */
  real fi, xi;
  int i;

  /* check that the center value is within bounds. */
  if (center < 0.0 || center > 1.0)
    throw("locus argument %.3f out of bounds [0,1]", center);

  /* check that the start value is within bounds. */
  if (start < 0.0 || start > 1.0)
    throw("start argument %.3f out of bounds [0,1]", start);

  /* check that the end value is within bounds. */
  if (end < 0.0 || end > 1.0)
    throw("end argument %.3f out of bounds [0,1]", end);

  /* ensure that the destination array is allocated. */
  if (!hx_window_alloc(wnd, d, len))
    throw("failed to allocate triangular array");

  /* compute the values of the window. */
  for (i = 0; i < len; i++) {
    /* compute the fractional index. */
    fi = ((real) i) / ((real) (len - 1));

    /* compute the current window value. */
    if (fi > center) {
      /* final: downward slope. */
      xi = fi * (end - 1.0) / (1.0 - center);
    }
    else {
      /* initial: upward slope. */
      xi = fi * (1.0 - start) / center;
    }

    /* store the computed window value. */
    wnd->x[i * wnd->n] = xi;
  }

  /* return success. */
  return 1;
}

/* hx_window_black(): allocate and construct a blackman window function.
 */
int hx_window_black (hx_array *wnd, int d, int len) {
  /* declare required variables:
   * @i: integer point index.
   * @fi: fractional point index.
   */
  real fi, xi;
  int i;

  /* ensure that the destination array is allocated. */
  if (!hx_window_alloc(wnd, d, len))
    throw("failed to allocate blackman array");

  /* compute the values of the window. */
  for (i = 0; i < len; i++) {
    /* compute the fractional index. */
    fi = ((real) i) / ((real) (len - 1));

    /* compute the current window value. */
    xi = 0.42 - 0.5 * cos(2.0 * M_PI * fi) + 0.08 * cos(4.0 * M_PI * fi);

    /* store the computed window value. */
    wnd->x[i * wnd->n] = xi;
  }

  /* return success. */
  return 1;
}

