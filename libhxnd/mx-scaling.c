
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

/* include the scaling method header. */
#include <hxnd/mx-scaling.h>

/* function declarations for all scaling methods:
 */
int mx_scale_none (hx_array *X);
int mx_scale_uv (hx_array *X);
int mx_scale_pareto (hx_array *X);
int mx_scale_range (hx_array *X);
int mx_scale_level (hx_array *X);
int mx_scale_vast (hx_array *X);

/* mx_scaling_func: function prototype for all scaling methods.
 * @X: pointer to the data matrix to scale.
 */
typedef int (*mx_scaling_func) (hx_array *X);

/* method_def: structure definition for looking up scaling method types by
 * string name.
 */
struct method_def {
  /* @name: scaling method string name.
   * @type: scaling method type.
   * @func: scaling function pointer.
   */
  const char *name;
  enum mx_scaling_type type;
  mx_scaling_func func;
};

/* methods: structure array of all supported scaling methods.
 */
const struct method_def methods[] = {
  { MX_SCALING_NAME_NONE,
    MX_SCALING_TYPE_NONE,
    &mx_scale_none
  },
  { MX_SCALING_NAME_UV,
    MX_SCALING_TYPE_UV,
    &mx_scale_uv
  },
  { MX_SCALING_NAME_PARETO,
    MX_SCALING_TYPE_PARETO,
    &mx_scale_pareto
  },
  { MX_SCALING_NAME_RANGE,
    MX_SCALING_TYPE_RANGE,
    &mx_scale_range
  },
  { MX_SCALING_NAME_LEVEL,
    MX_SCALING_TYPE_LEVEL,
    &mx_scale_level
  },
  { MX_SCALING_NAME_VAST,
    MX_SCALING_TYPE_VAST,
    &mx_scale_vast
  },

  /* null terminator. */
  { NULL, MX_SCALING_TYPE_UNDEFINED, NULL }
};

/* mx_scaling_lookup_type(): returns the enumerated type based on its string
 * representation.
 * @name: the scaling method name string.
 */
enum mx_scaling_type mx_scaling_lookup_type (const char *name) {
  /* declare a required variable. */
  unsigned int i;

  /* loop over all supported scaling methods. */
  for (i = 0; methods[i].name; i++) {
    /* check if the scaling method name matches. */
    if (strcmp(name, methods[i].name) == 0)
      return methods[i].type;
  }

  /* return an undefined method type. */
  return MX_SCALING_TYPE_UNDEFINED;
}

/* mx_scale(): apply a specified scaling method to a data matrix.
 * @X: pointer to the data matrix to scale.
 * @type: scaling method type.
 */
int mx_scale (hx_array *X, enum mx_scaling_type type) {
  /* declare a required variable. */
  unsigned int i;

  /* loop over all supported scaling methods. */
  for (i = 0; methods[i].name; i++) {
    /* check if the scaling method name matches. */
    if (methods[i].type == type)
      return methods[i].func(X);
  }

  /* return failure. */
  throw("unsupported scaling method");
}

/* mx_scale_by_name(): apply a specified (named) scaling method
 * to a data matrix.
 * @X: pointer to the data matrix to scale.
 * @name: scaling method string name.
 */
int mx_scale_by_name (hx_array *X, const char *name) {
  /* look up the specified scaling method and return its result. */
  return mx_scale(X, mx_scaling_lookup_type(name));
}

/* mx_scale_none(): mean-center the columns of a data matrix.
 * @X: data matrix to modify.
 */
int mx_scale_none (hx_array *X) {
  /* FIXME: implement mx_scale_none() */

  /* return success. */
  return 1;
}

/* mx_scale_uv(): mean-center and scale the columns of a data matrix
 * by their standard deviation.
 * @X: data matrix to modify.
 */
int mx_scale_uv (hx_array *X) {
  /* FIXME: implement mx_scale_uv() */

  /* return success. */
  return 1;
}

/* mx_scale_pareto(): mean-center and scale the columns of a data matrix
 * by the square root of their standard deviation.
 * @X: data matrix to modify.
 */
int mx_scale_pareto (hx_array *X) {
  /* FIXME: implement mx_scale_pareto() */

  /* return success. */
  return 1;
}

/* mx_scale_range(): mean-center and scale the columns of a data matrix
 * by their range.
 * @X: data matrix to modify.
 */
int mx_scale_range (hx_array *X) {
  /* FIXME: implement mx_scale_range() */

  /* return success. */
  return 1;
}

/* mx_scale_level(): mean-center and scale the columns of a data matrix
 * by their mean.
 * @X: data matrix to modify.
 */
int mx_scale_level (hx_array *X) {
  /* FIXME: implement mx_scale_level() */

  /* return success. */
  return 1;
}

/* mx_scale_vast(): mean-center and scale the columns of a data matrix
 * by their relative variance.
 * @X: data matrix to modify.
 */
int mx_scale_vast (hx_array *X) {
  /* FIXME: implement mx_scale_vast() */

  /* return success. */
  return 1;
}

