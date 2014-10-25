
/* ndmath: A framework for n-dimensional hypercomplex calculations for NMR.
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
#ifndef __NDMATH_HX_INDEX_H__
#define __NDMATH_HX_INDEX_H__

/* function declarations: */

int *hx_array_index_alloc (int k);

int *hx_array_index_build (int k, ...);

int hx_array_index_init (int *arr, int k);

int hx_array_index_pack (int k, int *sz, int *arr, int *pidx);

int hx_array_index_unpack (int k, int *sz, int **parr, int idx);

int hx_array_index_inc (int k, int *sz, int **parr);

int hx_array_index_diff (int k, int *a, int *b, int **pc);

int hx_array_index_bounded (int k, int *arr, int *lower, int *upper);

#endif /* __NDMATH_HX_INDEX_H__ */

