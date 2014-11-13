
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

/* fn_argdef_phase: define all accepted arguments for the 'phase' function.
 */
static const fn_args fn_argdef_phase[] = {
  { "ph0",     FN_ARGTYPE_FLOAT, "0.0" },
  { "ph1",     FN_ARGTYPE_FLOAT, "0.0" },
  { "inverse", FN_ARGTYPE_BOOL,  "0" },
  { NULL, '\0', NULL }
};

/* fn_execute_phase(): applies a single-dimension phase correction operation.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application.
 * @args: function argument string.
 */
int fn_execute_phase (datum *D, const int dim, const char *argstr) {
  /* declare variables to hold argument values.
   */
  real ph0, ph1;
  int inv;

  /* parse the function argument string. */
  if (!fn_scan_args(argstr, fn_argdef_phase, &ph0, &ph1, &inv))
    throw("failed to parse phase arguments");

  /* check the dimension index. */
  if (dim < 0 || dim >= D->nd)
    throw("dimension index %d out of bounds [0,%u)", dim, D->nd);

  /* FIXME: implement fn_execute_phase() */

  /* return success. */
  return 1;
}

