
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

/* fn_add(): adds a constant, or another file to a datum structure.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application, or -1.
 * @args: function argument definition array.
 */
int fn_add (datum *D, const int dim, const fn_arg *args) {
  /* declare variables to hold argument values.
   * @Dadd: datum structure (or NULL) of added hx-format file.
   * @fadd: filename string (or NULL) of added hx-format file.
   * @fscale: scaling factor for file addition.
   * @cadd: floating point real value to add.
   * @hxadd: hypercomplex scalar for constant addition.
   */
  real cadd, fscale;
  hx_scalar hxadd;
  int sub, d, n;
  datum Dadd;
  char *fadd;

  /* get the argument values from the argdef array. */
  if (!fn_args_get_all(args, &cadd, &fadd, &fscale, &sub))
    throw("failed to get add arguments");

  /* check if a subtraction was specified. */
  if (sub)
    fscale *= -1.0;

  /* check if a file argument was provided. */
  if (fadd) {
    /* attempt to load the file into a new datum. */
    if (!datum_load(&Dadd, fadd))
      throw("failed to read '%s'", fadd);

    /* add the two arrays together. */
    if (!hx_array_add_array(&D->array, &Dadd.array, fscale, &D->array))
      throw("failed to add file '%s'", fadd);

    /* free the temporary datum and the filename string. */
    datum_free(&Dadd);
    free(fadd);
  }

  /* check that the constant value is nonzero. */
  if (cadd) {
    /* store a local copy of the dimension index. */
    d = dim;

    /* check the dimension index (upper bound). */
    if (d >= D->nd)
      throw("dimension index %d out of bounds [0,%u)", d, D->nd);

    /* check the dimension index. */
    if (d < 0) {
      /* no dimension set: real addition. */
      d = 0;
      n = 0;
    }
    else {
      /* dimension set: complex addition. */
      n = 1 << d;
    }

    /* allocate a temporary scalar for the addition operation. */
    if (!hx_scalar_alloc(&hxadd, D->array.d))
      throw("failed to allocate scalar addition operand");

    /* set up the scalar value for addition. */
    hxadd.x[n] = cadd;

    /* perform the addition. */
    if (!hx_array_add_scalar(&D->array, &hxadd, 1.0, &D->array))
      throw("failed to add scalar value %f(%d)", cadd, n);

    /* free the allocated temporary scalar. */
    hx_scalar_free(&hxadd);
  }

  /* return success. */
  return 1;
}

