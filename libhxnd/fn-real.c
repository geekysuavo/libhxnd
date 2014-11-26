
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

/* fn_argdef_real: define all accepted arguments for the 'real' function.
 */
static const fn_args fn_argdef_real[] = {
  { NULL, '\0', NULL }
};

/* fn_execute_real(): drop imaginaries from the array of a datum structure.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application.
 * @args: function argument string.
 */
int fn_execute_real (datum *D, const int dim, const char *argstr) {
  /* base function behaviour on the dimension index. */
  if (dim < 0) {
    /* FIXME: implement complete real compression. */
  }
  else if (dim < D->nd) {
    /* FIXME: implement partial real compression. */
  }
  else
    throw("dimension index %d out of bounds [0,%u)", dim, D->nd);

  /* return success. */
  return 1;
}

