
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

/* ensure once-only inclusion. */
#ifndef __HXND_MX_DATASET_H__
#define __HXND_MX_DATASET_H__

/* include the nmr datum and function headers. */
#include <hxnd/nmr-datum.h>
#include <hxnd/fn.h>

/* dataset: data type for multiple data..
 *
 * dataset structures hold an array of datum structures and a corresponding
 * array of data handling functions that link the datum array to a data
 * matrix. thus, a dataset holds a dynamically updatable data matrix that
 * depends on (a) the input data and (b) the handling arguments.
 */
typedef struct {
  /* @D: array of datum structures.
   * @n: number of datum structures in the array..
   */
  datum *D;
  int n;

  /* @mask: array of integers indicating whether each datum structure
   * in the above array is to be included in the final data matrix.
   */
  int *mask;

  /* @funcs: list of dataset functions that transform the datum array
   * above (after masking) into the data matrix below.
   */
  fn_list funcs;

  /* @X: final real data matrix.
   * @X_ok: whether the data matrix is fresh.
   * @N: final data matrix observation (row) count.
   * @K: final data matrix variable (column) count.
   */
  int N, K, X_ok;
  hx_array X;
}
dataset;

/* function declarations (mx-dataset.c): */

void dataset_init (dataset *Dset);

void dataset_free (dataset *Dset);

int dataset_append_datum (dataset *Dset, datum *D);

int dataset_remove_datum (dataset *Dset, int idx);

void dataset_mask_datum (dataset *Dset, int idx);

void dataset_unmask_datum (dataset *Dset, int idx);

#endif /* __HXND_MX_DATASET_H__ */

