
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

/* include required standard c library headers. */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

/* hx_index: defined type holding an array of unpacked array indices.
 */
typedef int *hx_index;

/* define a print macro that writes the variable name. */
#define hx_index_print(k, arr) \
  hx_index_printfn(k, arr, #arr)

/* function declarations: */

hx_index hx_index_alloc (int k);

hx_index hx_index_build (int k, ...);

void hx_index_free (hx_index idx);

void hx_index_init (int k, hx_index idx);

hx_index hx_index_copy (int k, hx_index idx);

void hx_index_pack (int k, hx_index sz, hx_index idx, int *pidx);

void hx_index_unpack (int k, hx_index sz, hx_index idx, int pidx);

void hx_index_pack_tiled (int k, hx_index ntile, hx_index sztile,
                          hx_index idx, hx_index idxt, int *pidx);

int hx_index_incr (int k, hx_index sz, hx_index idx);

int hx_index_decr (int k, hx_index sz, hx_index idx);

int hx_index_incr_rev (int k, hx_index sz, hx_index idx);

int hx_index_decr_rev (int k, hx_index sz, hx_index idx);

int hx_index_incr_mask (int k, hx_index sz, hx_index idx, hx_index mask);

int hx_index_skip (int k, hx_index sz, hx_index idx, int kskip);

int hx_index_jump_init (int k, hx_index sz, int kskip,
                        int *ja, int *jb, int *jmax);

int hx_index_jump (int j, int ja, int jb);

int hx_index_diff (int k, hx_index a, hx_index b, hx_index c);

int hx_index_cmp (int k, hx_index a, hx_index b);

int hx_index_bounded (int k, hx_index idx, hx_index lower, hx_index upper);

int hx_index_sort (int k, hx_index idx);

hx_index hx_index_scheduled (int k, hx_index sz, int dsched, int nsched,
                             hx_index sched);

hx_index hx_index_unscheduled (int k, hx_index sz, int dsched, int nsched,
                               hx_index sched);

void hx_index_printfn (int k, hx_index idx, const char *s);

#endif /* __HXND_HX_INDEX_H__ */

