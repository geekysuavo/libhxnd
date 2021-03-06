
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

/* fn_subsamp(): nonuniformly subsamples a uniformly sampled datum.
 * @D: pointer to the datum to manipulate (in-place).
 * @dim: dimension of function application, or -1.
 * @args: function argument definition array.
 */
int fn_subsamp (datum *D, const int dim, const fn_arg *args) {
  /* declare variables to hold argument values.
   * @fsched: filename string of the input schedule.
   */
  char *fsched;

  /* declare a few required variables:
   * @idx: sampling schedule multidimensional index.
   * @pidx: packed linear sampling schedule index.
   * @d: datum dimension index.
   */
  hx_index idx, sz, zeros;
  int i, d, pidx, nzeros;

  /* get the argument values from the argdef array. */
  if (!fn_args_get_all(args, &fsched))
    throw("failed to get subsample arguments");

  /* check that a file argument was provided. */
  if (!fsched)
    throw("input schedule filename required");

  /* free the current datum schedule information. */
  if (!datum_sched_free(D))
    throw("failed to free current sampling schedule");

  /* load the new datum schedule information. */
  if (!datum_sched_read(D, fsched))
    throw("failed to read new sampling schedule from '%s'", fsched);

  /* ensure the loaded schedule has the correct dimensionality. */
  if (D->d_sched != D->nd - 1)
    throw("sampling schedule in '%s' has invalid dimensionality", fsched);

  /* allocate the schedule size array. */
  sz = hx_index_alloc(D->d_sched);

  /* check that allocation was successful. */
  if (!sz)
    throw("failed to allocate %d indices", D->d_sched);

  /* build the schedule size array. */
  for (i = 0, nzeros = 1; i < D->d_sched; i++) {
    /* store the current schedule size and increase the total count. */
    sz[i] = D->array.sz[i + 1];
    nzeros *= sz[i];
  }

  /* compute the number of unscheduled indices. */
  nzeros -= D->n_sched;

  /* create a list of unscheduled array indices. */
  zeros = hx_index_unscheduled(D->d_sched, sz, D->d_sched, D->n_sched,
                               D->sched);

  /* check that creation was successful. */
  if (!zeros)
    throw("failed to compute unscheduled array indices");

  /* allocate a multidimensional index. */
  idx = hx_index_alloc(D->nd);

  /* check that allocation was successful. */
  if (!idx)
    throw("failed to allocate %d indices", D->nd);

  /* loop over the points in the sampling schedule. */
  for (i = 0; i < nzeros; i++) {
    /* zero the index elements. */
    hx_index_init(D->nd, idx);

    /* unpack the currently unscheduled sub-index, then immediately
     * repack the linear index into a complete-array multidimensional
     * index.
     */
    hx_index_unpack(D->d_sched, sz, idx + 1, zeros[i]);
    hx_index_pack(D->array.k, D->array.sz, idx, &pidx);

    /* check that the linear index is in bounds. */
    if (pidx >= D->array.len / D->array.n)
      throw("sampling schedule entry #%d out of bounds", i);

    /* zero out the trace corresponding to the current index. */
    if (!hx_data_zero(D->array.x + pidx * D->array.n,
                      D->array.sz[0] * D->array.n))
      throw("failed to zero currently indexed trace");
  }

  /* set all indirect dimensions to nonuniformly sampled. */
  for (d = 1; d < D->nd; d++)
    D->dims[d].nus = 1;

  /* free the array indices and the filename string. */
  hx_index_free(zeros);
  hx_index_free(idx);
  hx_index_free(sz);
  free(fsched);

  /* return success. */
  return 1;
}

