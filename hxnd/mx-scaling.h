
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
#ifndef __HXND_MX_SCALING_H__
#define __HXND_MX_SCALING_H__

/* include the multivariate math header. */
#include <hxnd/mx.h>

/* define string constants for supported scaling method types.
 */
#define MX_SCALING_NAME_NONE    "none"
#define MX_SCALING_NAME_UV      "uv"
#define MX_SCALING_NAME_PARETO  "pareto"
#define MX_SCALING_NAME_RANGE   "range"
#define MX_SCALING_NAME_LEVEL   "level"
#define MX_SCALING_NAME_VAST    "vast"

/* mx_scaling_type: enumerated type for scaling method types.
 */
enum mx_scaling_type {
  MX_SCALING_TYPE_UNDEFINED,
  MX_SCALING_TYPE_NONE,
  MX_SCALING_TYPE_UV,
  MX_SCALING_TYPE_PARETO,
  MX_SCALING_TYPE_RANGE,
  MX_SCALING_TYPE_LEVEL,
  MX_SCALING_TYPE_VAST
};

/* function declarations: */

enum mx_scaling_type mx_scaling_lookup_type (const char *name);

int mx_scale (hx_array *X, enum mx_scaling_type type);

int mx_scale_by_name (hx_array *X, const char *name);

#endif /* __HXND_MX_SCALING_H__ */

