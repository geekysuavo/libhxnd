
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

/* fn_argdef_ht: define all accepted arguments for the 'ht' function.
 */
static const fn_args fn_argdef_ht[] = {
  { "zerofill", FN_ARGTYPE_BOOL, "0" },
  { "halve",    FN_ARGTYPE_BOOL, "0" },
  { "mirror",   FN_ARGTYPE_BOOL, "0" },
  { NULL, '\0', NULL }
};

/* fn_execute_ht(): applies a Hilbert transform to reconstruct imaginaries.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application.
 * @args: function argument string.
 */
int fn_execute_ht (datum *D, const int dim, const char *argstr) {
  /* declare variables to hold argument values.
   */
  int alt, neg;

  /* parse the function argument string. */
  if (!fn_scan_args(argstr, fn_argdef_ht, &alt, &neg))
    throw("failed to parse ht arguments");

  /* check the dimension index. */
  if (dim < 0 || dim >= D->nd)
    throw("dimension index %d out of bounds [0,%u)", dim, D->nd);

  /* FIXME: implement fn_execute_ht() */

  /* return success. */
  return 1;
}

