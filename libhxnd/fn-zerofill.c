
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

/* fn_argdef_zerofill: define all accepted arguments for the 'zerofill'
 * function.
 */
static const fn_args fn_argdef_zerofill[] = {
  { "times", FN_ARGTYPE_INT, "0" },
  { NULL, '\0', NULL }
};

/* fn_execute_zerofill(): zero-fill the array of a datum structure.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application.
 * @args: function argument string.
 */
int fn_execute_zerofill (datum *D, const int dim, const char *argstr) {
  /* declare variables to hold argument values.
   * @sznew: new size array for the core hypercomplex array.
   * @szv: new size values (with array count) for reshapes.
   * @szd: new size value for single-dimension operations.
   * @i: loop counter.
   * @k: topological dimension index.
   */
  int *sznew, nzf, i, k;
  unsigned int n, nx;

  /* parse the function argument string. */
  if (!fn_scan_args(argstr, fn_argdef_zerofill, &nzf))
    throw("failed to parse zerofill arguments");

  /* allocate the size array. */
  sznew = hx_array_index_alloc(D->array.k);;

  /* check that allocation succeeded. */
  if (!sznew)
    throw("failed to allocate new size array");

  /* compute the zero-filling multiplier. */
  for (i = 0, nx = 1; i < nzf; i++)
    nx *= 2;

  /* base function behaviour on the dimension index. */
  if (dim < 0) {
    /* zero-fill every dimension of the datum. */
    for (i = 0; i < D->array.k; i++) {
      /* compute the zero-filled value. */
      n = (unsigned int) D->array.sz[i];
      n = (hx_ispow2(n) ? n : hx_nextpow2(n));

      /* store the new value. */
      sznew[i] = n * nx;
    }
  }
  else if (dim < D->nd) {
    /* get the topological dimension index. */
    k = D->dims[dim].k;

    /* check that the index is in bounds. */
    if (k < 0 || k >= D->array.k)
      throw("topological dimension %d out of bounds [0,%d)", k, D->array.k);

    /* copy the sizes from the current array. */
    for (i = 0; i < D->array.k; i++)
      sznew[i] = D->array.sz[i];

    /* compute the zero-filled value of the indexed dimension. */
    n = sznew[k];
    n = (hx_ispow2(n) ? n : hx_nextpow2(n));

    /* set the new value. */
    sznew[k] = n * nx;
  }
  else
    throw("dimension index %d out of bounds [0,%u)", dim, D->nd);

  /* resize the datum content. */
  if (!datum_resize_array(D, sznew))
    throw("failed to resize datum");

  /* free the size array. */
  free(sznew);

  /* return success. */
  return 1;
}

