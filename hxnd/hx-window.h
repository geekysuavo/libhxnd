
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
#ifndef __HXND_HX_WINDOW_H__
#define __HXND_HX_WINDOW_H__

/* include the n-dimensional math header. */
#include <hxnd/hx.h>

/* define string constants for supported window function types.
 */
#define HX_WINDOW_NAME_SINE   "sine"
#define HX_WINDOW_NAME_EXP    "exp"
#define HX_WINDOW_NAME_GAUSS  "gauss"
#define HX_WINDOW_NAME_TRAP   "trap"
#define HX_WINDOW_NAME_TRI    "tri"
#define HX_WINDOW_NAME_BLACK  "black"

/* hx_window_type: enumerated type for apodization window types.
 */
enum hx_window_type {
  HX_WINDOW_TYPE_UNDEFINED,
  HX_WINDOW_TYPE_SINE,
  HX_WINDOW_TYPE_EXP,
  HX_WINDOW_TYPE_GAUSS,
  HX_WINDOW_TYPE_TRAP,
  HX_WINDOW_TYPE_TRI,
  HX_WINDOW_TYPE_BLACK
};

/* function declarations: */

enum hx_window_type hx_window_lookup_type (const char *name);

int hx_window_alloc (hx_array *wnd, int d, int len);

int hx_window_sine (hx_array *wnd, int d, int len, real width,
                    real start, real end, real order);

int hx_window_exp (hx_array *wnd, int d, int len, real width,
                   real lb);

int hx_window_gauss (hx_array *wnd, int d, int len, real width,
                     real invlb, real lb, real center);

int hx_window_trap (hx_array *wnd, int d, int len, real width,
                    real start, real end);

int hx_window_tri (hx_array *wnd, int d, int len, real width,
                   real center, real start, real end);

int hx_window_black (hx_array *wnd, int d, int len);

#endif /* __HXND_HX_WINDOW_H__ */

