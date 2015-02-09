
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
#ifndef __HXND_MX_STATS_H__
#define __HXND_MX_STATS_H__

/* include the n-dimensional math header. */
#include <hxnd/hx.h>

/* function declarations: */

int mx_stats_min (hx_array *x, int k, hx_array *m);

int mx_stats_max (hx_array *x, int k, hx_array *m);

int mx_stats_range (hx_array *x, int k, hx_array *r);

int mx_stats_med (hx_array *x, int k, hx_array *m);

int mx_stats_mean (hx_array *x, int k, hx_array *m);

int mx_stats_var (hx_array *x, int k, hx_array *s);

int mx_stats_stdev (hx_array *x, int k, hx_array *s);

#endif /* __HXND_MX_STATS_H__ */

