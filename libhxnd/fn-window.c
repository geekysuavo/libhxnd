
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

/* define all accepted window types.
 */
#define HX_WINTYPE_SINE
#define HX_WINTYPE_EXP
#define HX_WINTYPE_GAUSS
#define HX_WINTYPE_L2G
#define HX_WINTYPE_TRAP

/* fn_argdef_window: define all accepted arguments for the 'window' function.
 */
static const fn_args fn_argdef_window[] = {
  { "type", FN_ARGTYPE_STRING, NULL },
  /* FIXME: define all 'window' arguments. */
  { NULL, '\0', NULL }
};

/* fn_execute_widow(): applies a window function to a datum structure.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application.
 * @args: function argument string.
 */
int fn_execute_window (datum *D, const int dim, const char *argstr) {
  /* declare variables to hold argument values.
   */
  char *wtype;

  /* parse the function argument string. */
  if (!fn_scan_args(argstr, fn_argdef_window, &wtype))
    throw("failed to parse window arguments");

  /* FIXME: implement the 'window' function. */

  /* return success. */
  return 1;
}

