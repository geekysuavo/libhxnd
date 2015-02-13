
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

/* include the processing function header. */
#include <hxnd/fn.h>
#include <hxnd/fn-handlers.h>

/* define type strings for all available projection operations. */
#define FN_PROJECT_SUM  "sum"
#define FN_PROJECT_MIN  "min"
#define FN_PROJECT_MAX  "max"

/* fn_project_sum(): projector callback for sum-type projection.
 * see fn_project() and hx_array_projector_cb() for more details.
 */
int fn_project_sum (hx_array *y, real *val) {
  /* declare a required index variable. */
  int i;

  /* initialize the projection sum. */
  hx_data_zero(val, y->n);

  /* loop over the vector elements. */
  for (i = 0; i < y->len; i += y->n)
    hx_data_add(y->x + i, val, val, 1.0, y->d, y->n);

  /* return success. */
  return 1;
}

/* fn_project_min(): projector callback for floor-type projection.
 * see fn_project() and hx_array_projector_cb() for more details.
 */
int fn_project_min (hx_array *y, real *val) {
  /* declare a few required variables. */
  real yi, ymax;
  int i, imax;

  /* loop over the vector elements. */
  for (i = imax = 0, ymax = 0.0; i < y->len; i += y->n) {
    /* compute the current vector element norm. */
    yi = hx_data_real_norm(y->x + i, y->d, y->n);

    /* check if the current norm is greater than previously found. */
    if (yi > ymax) {
      /* store the new index. */
      imax = i;
      ymax = yi;
    }
  }

  /* copy the data from the identified maximum element. */
  hx_data_copy(y->x + imax, val, y->n);

  /* return success. */
  return 1;
}

/* fn_project_max(): projector callback for skyline-type projection.
 * see fn_project() and hx_array_projector_cb() for more details.
 */
int fn_project_max (hx_array *y, real *val) {
  /* declare a few required variables. */
  real yi, ymin;
  int i, imin;

  /* compute the first vector element norm. */
  ymin = hx_data_real_norm(y->x, y->d, y->n);
  imin = 0;

  /* loop over the vector elements. */
  for (i = 1; i < y->len; i += y->n) {
    /* compute the current vector element norm. */
    yi = hx_data_real_norm(y->x + i, y->d, y->n);

    /* check if the current norm is less than previously found. */
    if (yi < ymin) {
      /* store the new index. */
      imin = i;
      ymin = yi;
    }
  }

  /* copy the data from the identified minimum element. */
  hx_data_copy(y->x + imin, val, y->n);

  /* return success. */
  return 1;
}

/* fn_project(): project a dimension out of the array of a datum structure.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application, or -1.
 * @args: function argument definition array.
 */
int fn_project (datum *D, const int dim, const fn_arg *args) {
  /* declare variables to hold argument values. */
  char *type = NULL;

  /* declare a required variable:
   * @projector: projector callback function to use.
   */
  hx_array_projector_cb projector = &fn_project_sum;

  /* get the argument values from the argdef array. */
  if (!fn_args_get_all(args, &type))
    throw("failed to get project arguments");

  /* check the dimension index. */
  if (dim < 0 || dim >= D->nd)
    throw("dimension index %d out of bounds [0,%u)", dim, D->nd);

  /* check if a projection type was specified. */
  if (type) {
    /* determine which projection type was specified. */
    if (strcmp(type, FN_PROJECT_SUM) == 0) {
      /* set the 'sum' function pointer. */
      projector = &fn_project_sum;
    }
    else if (strcmp(type, FN_PROJECT_MIN) == 0) {
      /* set the 'min' function pointer. */
      projector = &fn_project_min;
    }
    else if (strcmp(type, FN_PROJECT_MAX) == 0) {
      /* set the 'max' function pointer. */
      projector = &fn_project_max;
    }
    else
      throw("unsupported projection type '%s'", type);
  }

  /* execute the datum projection function. */
  if (!datum_array_project(D, dim, projector))
    throw("failed to project datum array");

  /* free the string argument value. */
  if (type)
    free(type);

  /* return success. */
  return 1;
}

