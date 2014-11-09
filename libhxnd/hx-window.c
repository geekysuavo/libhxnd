
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
  { HX_WINDOW_NAME_SINE,  HX_WINDOW_TYPE_SINE },
  { HX_WINDOW_NAME_EXP,   HX_WINDOW_TYPE_EXP },
  { HX_WINDOW_NAME_GAUSS, HX_WINDOW_TYPE_GAUSS },
  { HX_WINDOW_NAME_TRAP,  HX_WINDOW_TYPE_TRAP },
  { HX_WINDOW_NAME_TRI,   HX_WINDOW_TYPE_TRI },
  { NULL,                 HX_WINDOW_TYPE_UNDEFINED }
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
 * @width: spectral width parameter.
 */
int hx_window_alloc (hx_array *wnd, int d, int len, real width) {
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
  /* ensure that the destination array is allocated. */
  if (!hx_window_alloc(wnd, d, len, width))
    throw("failed to allocate sinusoidal array");

  /* FIXME: implement hx_window_sine() */

  /* return success. */
  return 1;
}

/* hx_window_exp(): allocate and construct an exponential window function.
 * @lb: exponential line-broadening in Hertz.
 */
int hx_window_exp (hx_array *wnd, int d, int len, real width,
                   real lb) {
  /* ensure that the destination array is allocated. */
  if (!hx_window_alloc(wnd, d, len, width))
    throw("failed to allocate exponential array");

  /* FIXME: implement hx_window_exp() */

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
  /* ensure that the destination array is allocated. */
  if (!hx_window_alloc(wnd, d, len, width))
    throw("failed to allocate gaussian array");

  /* FIXME: implement hx_window_gauss() */

  /* return success. */
  return 1;
}

/* hx_window_trap(): allocate and construct a trapezoidal window function.
 * @start: size of starting ramp.
 * @end: size of ending ramp.
 */
int hx_window_trap (hx_array *wnd, int d, int len, real width,
                    real start, real end) {
  /* ensure that the destination array is allocated. */
  if (!hx_window_alloc(wnd, d, len, width))
    throw("failed to allocate trapezoidal array");

  /* FIXME: implement hx_window_trap() */

  /* return success. */
  return 1;
}

/* hx_window_tri(): allocate and construct a triangular window function.
 * @locus: position of triangle maximum.
 * @start: height of starting point.
 * @end: height of ending point.
 */
int hx_window_tri (hx_array *wnd, int d, int len, real width,
                   real locus, real start, real end) {
  /* ensure that the destination array is allocated. */
  if (!hx_window_alloc(wnd, d, len, width))
    throw("failed to allocate triangular array");

  /* FIXME: implement hx_window_tri() */

  /* return success. */
  return 1;
}

