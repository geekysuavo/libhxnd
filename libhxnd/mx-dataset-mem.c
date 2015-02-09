
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

/* include the nmr dataset header. */
#include <hxnd/mx-dataset.h>

/* dataset_init(): initializes the elements of a dataset structure.
 * @Dset: pointer to the dataset to initialize.
 */
void dataset_init (dataset *Dset) {
  /* initialize the datum array. */
  Dset->D = NULL;
  Dset->n = 0;

  /* initialize the mask. */
  Dset->mask = NULL;

  /* initialize the function list. */
  fn_list_init(&Dset->funcs);

  /* initialize the normalization and scaling arrays. */
  hx_array_init(&Dset->centers);
  hx_array_init(&Dset->scales);
  hx_array_init(&Dset->norms);

  /* indicate that the data matrix is not allocated. */
  Dset->X_ok = Dset->X_damage = 0;
  Dset->N = Dset->K = 0;
  hx_array_init(&Dset->X);
}

/* dataset_free(): frees all allocated innards of a dataset structure.
 * @Dset: pointer to the dataset to free.
 */
void dataset_free (dataset *Dset) {
  /* free the datum array, but do not free each datum. */
  if (Dset->D)
    free(Dset->D);

  /* free the mask. */
  if (Dset->mask)
    free(Dset->mask);

  /* free the function list. */
  fn_list_free(&Dset->funcs);

  /* free the data matrix. */
  dataset_matrix_free(Dset);

  /* re-initialize the dataset. */
  dataset_init(Dset);
}

