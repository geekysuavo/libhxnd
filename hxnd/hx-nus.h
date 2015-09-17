
/* hxnd: A framework for n-dimensional hypercomplex calculations for NMR.
 * Copyright (C) 2014-2015  Bradley Worley  <geekysuavo@gmail.com>.
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
#ifndef __HXND_HX_NUS_H__
#define __HXND_HX_NUS_H__

/* function declarations: */

int hx_array_ist (hx_array *x, hx_index dx, hx_index kx,
                  int dsched, int nsched, hx_index sched,
                  int niter, real thresh);

int hx_array_irls (hx_array *x, hx_index dx, hx_index kx,
                   int dsched, int nsched, hx_index sched,
                   int niter, real pa, real pb);

int hx_array_ffm (hx_array *x, hx_index dx, hx_index kx,
                  int dsched, int nsched, hx_index sched,
                  int niter, enum hx_entropy_type type);

#endif /* __HXND_HX_NUS_H__ */

