
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

/* hx_window_alloc(): allocates memory for a hypercomplex window function
 * and stores the computed values in the allocated array.
 * @wnd: hypercomplex array pointer for the result.
 * @d: dimension count for the array.
 * @len: length of the window vector.
 * @width: spectral width parameter.
 * @type: type of the window function.
 */
int hx_window_alloc (hx_array *wnd, int d, int len, real width,
                     enum hx_window_type type) {
  /* FIXME */

  /* return success. */
  return 1;
}

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

