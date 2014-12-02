
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
#ifndef __HXND_HX_INDEX_H__
#define __HXND_HX_INDEX_H__

/* define a print macro that writes the variable name. */
#define hx_array_index_print(k, arr) \
  hx_array_index_printfn(k, arr, #arr)

/* function declarations: */

int *hx_array_index_alloc (int k);

int *hx_array_index_build (int k, ...);

int hx_array_index_init (int *arr, int k);

int hx_array_index_pack (int k, int *sz, int *arr, int *pidx);

int hx_array_index_unpack (int k, int *sz, int *arr, int idx);

int hx_array_index_incr (int k, int *sz, int *arr);

int hx_array_index_decr (int k, int *sz, int *arr);

int hx_array_index_incr_rev (int k, int *sz, int *arr);

int hx_array_index_decr_rev (int k, int *sz, int *arr);

int hx_array_index_skip (int k, int *sz, int *arr, int kskip);

int hx_array_index_jump_init (int k, int *sz, int kskip,
                              int *ja, int *jb, int *jmax);

int hx_array_index_jump (int j, int ja, int jb);

int hx_array_index_diff (int k, int *a, int *b, int *c);

int hx_array_index_bounded (int k, int *arr, int *lower, int *upper);

int hx_array_index_sort (int k, int *order);

int *hx_array_index_scheduled (int k, int *sz, int dsched, int nsched,
                               int *sched);

int *hx_array_index_unscheduled (int k, int *sz, int dsched, int nsched,
                                 int *sched);

void hx_array_index_printfn (int k, int *arr, const char *s);

#endif /* __HXND_HX_INDEX_H__ */

