
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

/* fn_argdef_cut: define all accepted arguments for the 'cut' function.
 */
static const fn_args fn_argdef_cut[] = {
  { "trace", FN_ARGTYPE_INTS, NULL },
  { "plane", FN_ARGTYPE_INTS, NULL },
  { NULL, '\0', NULL }
};

/* fn_execute_cut(): cut a trace or a plane from the array of a
 * datum structure.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application.
 * @args: function argument string.
 */
int fn_execute_cut (datum *D, const int dim, const char *argstr) {
  /* declare variables to hold argument values.
   * @ivtr: int-array for extracting traces.
   * @ivpl: int-array for extracting planes.
   */
  int *ivtr, *ivpl;

  /* declare a few required variables:
   * @iv: int-array that points to either @ivtr or @ivpl.
   * @lower: lower bound index array for slicing.
   * @upper: upper bound index array for slicing.
   * @i: general-purpose loop counter.
   * @nz: number of zeros in @ivtr or @ivpl.
   */
  int *iv, *lower, *upper, i, nz;

  /* check that no dimension was specified. */
  if (dim >= 0)
    throw("dimension index specification not supported");

  /* parse the function argument string. */
  if (!fn_scan_args(argstr, fn_argdef_cut, &ivtr, &ivpl))
    throw("failed to parse cut arguments");

  /* ensure that both modes were not specified. */
  if (ivtr && ivpl) {
    /* free the arrays. */
    free(ivtr);
    free(ivpl);

    /* throw an exception. */
    throw("trace and plane cut modes are mutually exclusive");
  }

  /* ensure that at least one mode was specified. */
  if (!ivtr && !ivpl)
    throw("no cut mode specified");

  /* store the int-array in a variable that reduces the amount of
   * code duplication below.
   */
  iv = (ivtr ? ivtr : ivpl);

  /* allocate the lower and upper bound index arrays. */
  lower = hx_array_index_alloc(D->array.k);
  upper = hx_array_index_alloc(D->array.k);

  /* check that allocation succeeded. */
  if (!lower || !upper)
    throw("failed to allocate index arrays");

  /* check that the int-array is of correct length. */
  if (iv[0] != D->array.k)
    throw("invalid array length (%d != %d)", iv[0], D->array.k);

  /* build the lower and upper bound arrays. */
  for (i = 0, nz = 0; i < D->array.k; i++) {
    /* store the lower bound. */
    lower[i] = (iv[i + 1] > 0 ? iv[i + 1] - 1 : 0);

    /* store the upper bound. */
    upper[i] = (iv[i + 1] > 0 ? iv[i + 1] - 1 : D->array.sz[i] - 1);

    /* count the number of zero-valued indices. */
    if (iv[i + 1] < 1)
      nz++;
  }

  /* check that the index array contains only one invalid entry. */
  if (ivtr && nz != 1)
    throw("trace cutting requires one zero-valued index");
  else if (ivpl && nz != 2)
    throw("plane cutting requires two zero-valued indices");

  /* free the trace/plane int-array. */
  free(iv);

  /* check all values in the built lower and upper bounds. */
  for (i = 0; i < D->array.k; i++) {
    /* check the lower bound. */
    if (lower[i] < 0 || lower[i] >= D->array.sz[i]) {
      throw("lower bound %d (#%d) out of range [0,%d)",
            lower[i], i, D->array.sz[i]);
    }

    /* check the upper bound. */
    if (upper[i] < 0 || upper[i] >= D->array.sz[i]) {
      throw("upper bound %d (#%d) out of range [0,%d)",
            upper[i], i, D->array.sz[i]);
    }
  }

  /* slice the datum based on the built lower and upper bounds. */
  if (!datum_array_slice(D, lower, upper))
    throw("failed to perform slice operation");

  /* free the lower and upper bound index arrays. */
  free(lower);
  free(upper);

  /* return success. */
  return 1;
}

