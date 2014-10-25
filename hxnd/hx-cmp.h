
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
#ifndef __HXND_HX_CMP_H__
#define __HXND_HX_CMP_H__

/* definitions of comparison mismatch codes.
 */
#define HXCMP_ID    0  /* identity, equality. */
#define HXCMP_DIMS  1  /* algebraic dimensionality mismatch. */
#define HXCMP_TOPO  2  /* topological dimensionality mismatch. */
#define HXCMP_SIZE  3  /* topological size mismatch. */
#define HXCMP_DATA  4  /* coefficient data mismatch. */

/* function declarations, scalars: */

int hx_scalar_cmp (hx_scalar *a, hx_scalar *b);

int hx_scalar_dims_cmp (hx_scalar *a, hx_scalar *b);

/* function declarations, arrays: */

int hx_array_cmp (hx_array *a, hx_array *b);

int hx_array_dims_cmp (hx_array *a, hx_array *b);

int hx_array_topo_cmp (hx_array *a, hx_array *b);

int hx_array_conf_cmp (hx_array *a, hx_array *b);

#endif /* __HXND_HX_CMP_H__ */

