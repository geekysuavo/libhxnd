
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

/* include the window function header. */
#include <hxnd/hx-window.h>

/* fn_window(): applies a window function to a datum structure.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application, or -1.
 * @args: function argument definition array.
 */
int fn_window (datum *D, const int dim, const fn_arg *args) {
  /* declare a few required variables:
   * @type: the window function enumerated type to apply.
   * @ldim: a local (thus mutable) copy of the @dim value.
   * @d: algebraic dimensionality of the datum array.
   * @n: number of characters in the window type string.
   * @len: size of the windowed array dimension.
   * @ret: return value of the window construction.
   * @wnd: vector-shaped array of window coefficients.
   * @width: spectral width of the windowed dimension.
   */
  enum hx_window_type type;
  int ldim, d, n, len, ret;
  hx_array wnd;
  real width;

  /* declare variables to hold argument values. */
  real start, end, order, lb, invlb, center;
  char *stype;

  /* get the argument values from the argdef array. */
  if (!fn_args_get_all(args, &stype, &start, &end, &order,
                       &lb, &invlb, &center))
    throw("failed to get window arguments");

  /* store the dimensionality into a local variable. */
  ldim = dim;
  if (ldim < 0)
    ldim = 0;

  /* check the dimension index. */
  if (ldim >= D->nd)
    throw("dimension index %d out of bounds [0,%d)", ldim, D->nd);

  /* check that a type was provided. */
  if (!stype) {
    /* none provided. set the default window type. */
    n = strlen(HX_WINDOW_NAME_SINE) + 1;
    stype = (char*) malloc(n * sizeof(char));

    /* check that allocation was successful. */
    if (!stype)
      throw("failed to allocate window type string");

    /* store the default window name. */
    strcpy(stype, HX_WINDOW_NAME_SINE);
  }

  /* determine the window enumerated type. */
  type = hx_window_lookup_type(stype);

  /* store local window variables. */
  width = D->dims[ldim].width;
  len = D->array.sz[ldim];
  d = D->array.d;

  /* initialize the window array contents. */
  hx_array_init(&wnd);

  /* determine which window function to construct. */
  ret = 0;
  switch (type) {
    /* sine. */
    case HX_WINDOW_TYPE_SINE:
      /* construct a sine window. */
      ret = hx_window_sine(&wnd, d, len, width, start, end, order);
      break;

    /* exponential. */
    case HX_WINDOW_TYPE_EXP:
      /* construct an exponential window. */
      ret = hx_window_exp(&wnd, d, len, width, lb);
      break;

    /* gaussian. */
    case HX_WINDOW_TYPE_GAUSS:
      /* construct a gaussian window. */
      ret = hx_window_gauss(&wnd, d, len, width, invlb, lb, center);
      break;

    /* trapezoidal. */
    case HX_WINDOW_TYPE_TRAP:
      /* construct a trapezoidal window. */
      ret = hx_window_trap(&wnd, d, len, width, start, end);
      break;

    /* triangular. */
    case HX_WINDOW_TYPE_TRI:
      /* construct a triangular window. */
      ret = hx_window_tri(&wnd, d, len, width, center, start, end);
      break;

    /* blackman. */
    case HX_WINDOW_TYPE_BLACK:
      ret = hx_window_black(&wnd, d, len);
      break;

    /* undefined. */
    default:
    case HX_WINDOW_TYPE_UNDEFINED:
      throw("window type '%s' undefined", stype);
  }

  /* check if the array was successfully constructed. */
  if (!ret)
    throw("failed to construct %s window", stype);

  /* perform a trace-wise multiplication by the window. */
  if (!hx_array_mul_vector(&D->array, &wnd, D->dims[ldim].k, &D->array))
    throw("failed to perform window multiplication");

  /* free the allocated array. */
  hx_array_free(&wnd);

  /* free the window type string. */
  free(stype);

  /* return success. */
  return 1;
}

