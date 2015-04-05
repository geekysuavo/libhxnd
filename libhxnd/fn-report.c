
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

/* fn_report_sumsq(): report the sum of squared elements of a datum array.
 * @D: pointer to the datum to access.
 */
int fn_report_sumsq (datum *D) {
  /* declare a few required variables:
   * @pidx: packed linear array index.
   * @S: total sum of squares of the data.
   */
  int pidx;
  real S;

  /* loop over the datum array. */
  for (pidx = 0, S = 0.0; pidx < D->array.len; pidx += D->array.n)
    S += hx_data_real_norm(D->array.x + pidx, D->array.d, D->array.n);

  /* print the output value. */
  fprintf(stdout, "sumsq = %18.8le\n", S);

  /* return success. */
  return 1;
}

/* fn_report(): report specified statistics pertaining to the array of
 * a datum structure.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application, or -1.
 * @args: function argument definition array.
 */
int fn_report (datum *D, const int dim, const fn_arg *args) {
  /* declare variables to hold argument values. */
  int sumsq;

  /* get the argument values from the argdef array */
  if (!fn_args_get_all(args, &sumsq))
    throw("failed to get report arguments");

  /* check that no dimension was specified. */
  if (dim >= 0)
    throw("dimension index specification not supported");

  /* check if the sum of squares was requested. */
  if (sumsq && !fn_report_sumsq(D))
    throw("failed to report sum of squares");

  /* return success. */
  return 1;
}

