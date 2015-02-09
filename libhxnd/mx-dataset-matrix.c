
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

/* dataset_matrix_damage(): indicate that the data matrix contained within
 * a dataset structure has been damaged and requires rebuilding.
 * @Dset: pointer to the dataset to update.
 * @idmg: damage index to report.
 */
int dataset_matrix_damage (dataset *Dset, int idmg) {
  /* check that the damage index is in bounds. */
  if (idmg > Dset->funcs.n)
    return 1;

  /* report the damage index and flag the matrix as damaged. */
  Dset->X_damage = idmg;
  Dset->X_ok = 0;

  /* return success. */
  return 1;
}

/* dataset_matrix_build(): construct a data matrix from an array of
 * datum structures based on a defined set of functions. if the data
 * matrix already exists, destroy it and rebuild it from scratch.
 * @Dset: pointer to the dataset to update.
 */
int dataset_matrix_build (dataset *Dset) {
  /* declare a few required variables:
   * @ifn: function list index.
   */
  int ifn;

  /* return if the matrix is undamaged. */
  if (Dset->X_ok)
    return 1;

  /* loop from the point of damage to the function list end. */
  for (ifn = Dset->X_damage; ifn < Dset->funcs.n; ifn++) {
    /* execute the currently indexed dataset function. */
    if (!fn_execute(Dset, 0, &Dset->funcs.v[ifn], Dset->funcs.v[ifn].args))
      throw("failed to apply function '%s'", Dset->funcs.v[ifn].name);
  }

  /* flag the matrix as undamaged. */
  Dset->X_damage = Dset->funcs.n;
  Dset->X_ok = 1;

  /* return success. */
  return 1;
}

/* dataset_matrix_free(): completely free the allocated data matrix of
 * a dataset structure, including all related scaling and normalization
 * information.
 * @Dset: pointer to the dataset to modify.
 */
int dataset_matrix_free (dataset *Dset) {
  /* free the scaling and normalization arrays. */
  hx_array_free(&Dset->centers);
  hx_array_free(&Dset->scales);
  hx_array_free(&Dset->norms);

  /* free the data matrix. */
  hx_array_free(&Dset->X);

  /* report the data matrix as damaged. */
  dataset_matrix_damage(Dset, 0);

  /* return success. */
  return 1;
}

