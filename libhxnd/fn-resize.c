
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

/* fn_argdef_resize: define all accepted arguments for the 'resize' function.
 */
static const fn_args fn_argdef_resize[] = {
  { "size",  FN_ARGTYPE_INT,  "0" },
  { "shape", FN_ARGTYPE_INTS, NULL },
  { NULL, '\0', NULL }
};

/* fn_execute_resize(): resize or reshape the array of a datum structure.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application.
 * @args: function argument string.
 */
int fn_execute_resize (datum *D, const int dim, const char *argstr) {
  /* declare variables to hold argument values.
   * @sznew: new size array for the core hypercomplex array.
   * @szv: new size values (with array count) for reshapes.
   * @szd: new size value for single-dimension operations.
   * @i: loop counter.
   */
  int *sznew, *szv, szd, i;

  /* parse the function argument string. */
  if (!fn_scan_args(argstr, fn_argdef_resize, &szd, &szv))
    throw("failed to parse resize arguments");

  /* initialize the size array. */
  sznew = NULL;

  /* base function behaviour on the dimension index. */
  if (dim < 0) {
    /* check that the 'reshape' parameter was passed. */
    if (!szv)
      throw("expected 'shape' argument not found");

    /* check that the parsed int-array is of correct length. */
    if (szv[0] != D->array.k)
      throw("invalid array length (%d != %d)", szv[0], D->array.k);

    /* allocate memory for the size array. */
    sznew = hx_array_index_alloc(D->array.k);

    /* check that allocation succeeded. */
    if (!sznew)
      throw("failed to allocate new size array");

    /* copy the sizes from the parsed int-array. */
    for (i = 0; i < D->array.k; i++)
      sznew[i] = szv[i + 1];

    /* free the parsed int-array. */
    free(szv);
  }
  else if (dim < D->nd) {
    /* check that the 'resize' parameter was passed. */
    if (!szd)
      throw("expected 'size' argument not found");

    /* allocate memory for the size array. */
    sznew = hx_array_index_alloc(D->array.k);

    /* check that allocation succeeded. */
    if (!sznew)
      throw("failed to allocate new size array");

    /* copy the sizes from the current array. */
    for (i = 0; i < D->array.k; i++)
      sznew[i] = D->array.sz[i];

    /* get the topological dimension index. */
    i = D->dims[dim].k;

    /* check that the index is in bounds. */
    if (i < 0 || i >= D->array.k)
      throw("topological index %d out of bounds [0,%d)", i, D->array.k);

    /* set the modified size. */
    sznew[i] = szd;
  }
  else
    throw("dimension index %d out of bounds [0,%u)", dim, D->nd);

  /* resize the datum content. */
  if (!datum_array_resize(D, sznew))
    throw("failed to resize datum");

  /* free the size array. */
  free(sznew);

  /* return success. */
  return 1;
}

