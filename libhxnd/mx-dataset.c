
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

/* dataset_append_datum(): copy the contents of a datum structure into
 * the datum array of a dataset.
 * @Dset: pointer to the dataset to update.
 * @D: pointer to the datum to append.
 */
int dataset_append_datum (dataset *Dset, datum *D) {
  /* FIXME: implement dataset_append_datum() */
  throw("UNIMPLEMENTED");

  /* flag the data matrix as damaged. */
  dataset_matrix_damage(Dset, 0);

  /* return success. */
  return 1;
}

/* dataset_remove_datum(): remove an indexed datum from the datum
 * array of a dataset.
 * @Dset: pointer to the dataset to update.
 * @idx: array index of the datum to remove.
 */
int dataset_remove_datum (dataset *Dset, int idx) {
  /* check that the datum index is in bounds. */
  if (idx < 0 || idx >= Dset->n)
    throw("datum index %d out of bounds [0,%d)", idx, Dset->n);

  /* FIXME: implement dataset_remove_datum() */
  throw("UNIMPLEMENTED");

  /* flag the data matrix as damaged. */
  dataset_matrix_damage(Dset, 0);

  /* return success. */
  return 1;
}

/* dataset_mask_datum(): mask an indexed datum in a dataset.
 * @Dset: pointer to the dataset to update.
 * @idx: array index of the datum to mask.
 */
int dataset_mask_datum (dataset *Dset, int idx) {
  /* check that the datum index is in bounds. */
  if (idx < 0 || idx >= Dset->n)
    throw("datum index %d out of bounds [0,%d)", idx, Dset->n);

  /* return if the datum is already masked. */
  if (!Dset->mask[idx])
    return 1;

  /* disable the indexed datum in the mask array and return. */
  Dset->mask[idx] = 0;

  /* flag the data matrix as damaged. */
  dataset_matrix_damage(Dset, 0);

  /* return success. */
  return 1;
}

/* dataset_unmask_datum(): unmask an indexed datum in a dataset.
 * @Dset: pointer to the dataset to update.
 * @idx: array index of the datum to mask.
 */
int dataset_unmask_datum (dataset *Dset, int idx) {
  /* check that the datum index is in bounds. */
  if (idx < 0 || idx >= Dset->n)
    throw("datum index %d out of bounds [0,%d)", idx, Dset->n);

  /* return if the datum is already unmasked. */
  if (Dset->mask[idx])
    return 1;

  /* enable the indexed datum in the mask array. */
  Dset->mask[idx] = 1;

  /* flag the data matrix as damaged. */
  dataset_matrix_damage(Dset, 0);

  /* return success. */
  return 1;
}

