
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

/* fn_argdef_ist: define all accepted arguments for the 'ist' function.
 */
static const fn_args fn_argdef_ist[] = {
  { "thresh", FN_ARGTYPE_FLOAT, "0.9" },
  { "iters",  FN_ARGTYPE_INT,   "200" },
  { NULL, '\0', NULL }
};

/* fn_execute_ist(): reconstructs a nonuniformly sampled dimension using
 * iterative soft thresholding.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application.
 * @args: function argument string.
 */
int fn_execute_ist (datum *D, const int dim, const char *argstr) {
  /* declare variables to hold argument values. */
  real thresh;
  int iters;

  /* declare a few required variables. */
  unsigned int d;
  int nnus;

  /* parse the function argument string. */
  if (!fn_scan_args(argstr, fn_argdef_ist, &thresh, &iters))
    throw("failed to parse ist arguments");

  /* check that no dimension was specified. */
  if (dim >= 0)
    throw("dimension index specification not supported");

  /* count the number of nonuniform dimensions. */
  for (d = 0, nnus = 0; d < D->nd; d++)
    nnus += (D->dims[d].nus ? 1 : 0);

  /* check that the nonuniform dimension indices match expectations:
   * a. first dimension is uniformly sampled.
   * b. all other dimensions are subsampled.
   */
  if (D->dims[0].nus || nnus != D->nd - 1)
    throw("unexpected initial conditions for nus reconstruction");

  /* check that the nonuniform dimension count matches the schedule. */
  if (nnus != D->d_sched)
    throw("unexpected nus dimension count (%d != %d)", nnus, D->d_sched);

  /* execute the reconstruction. */
  if (!hx_array_ist(&D->array, D->d_sched, D->n_sched,
                    D->sched, iters, thresh))
    throw("failed to perform ist reconstruction");

  /* return success. */
  return 1;
}

