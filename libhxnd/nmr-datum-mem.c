
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

/* include the nmr data header. */
#include <hxnd/nmr-datum.h>

/* datum_init(): initializes the elements of an NMR datum structure.
 * @D: pointer to the datum to initialize.
 */
void datum_init (datum *D) {
  /* initialize the input data filename. */
  D->fname = NULL;

  /* initialize the datum type. */
  D->endian = BYTES_ENDIAN_AUTO;
  D->type = DATUM_TYPE_UNDEFINED;

  /* initialize the datum date. */
  D->epoch = time(NULL);

  /* initialize the dimension count and dimensions array. */
  D->dims = NULL;
  D->nd = 0;

  /* initialize the sampling schedule. */
  D->sched = NULL;
  D->d_sched = 0;
  D->n_sched = 0;

  /* initialize the group delay value. */
  D->grpdelay = 0.0;

  /* indicate that the array is not allocated. */
  D->array_alloc = 0;
  hx_array_init(&D->array);
}

/* datum_free(): frees all allocated innards of an NMR datum structure.
 * @D: pointer to the datum to free.
 */
void datum_free (datum *D) {
  /* free the input filename. */
  if (D->fname)
    free(D->fname);

  /* free the dimensions array. */
  if (D->dims)
    free(D->dims);

  /* free the sampling schedule. */
  if (D->sched)
    free(D->sched);

  /* free the array data. */
  datum_array_free(D);

  /* re-initialize the datum. */
  datum_init(D);
}

