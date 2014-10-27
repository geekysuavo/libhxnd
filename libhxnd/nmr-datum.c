
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

/* datum_print(): prints the metadata associated with an acquired NMR datum.
 * @D: the datum to print data from.
 * @fname: the output filename.
 */
int datum_print (datum *D, const char *fname) {
  /* declare a few required variables:
   * @d: dimension loop counter.
   * @fh: the file handle used for writing.
   */
  unsigned int d, n;
  FILE *fh;

  /* open the output file. */
  if (fname)
    fh = fopen(fname, "wb");
  else
    fh = stdout;

  /* check that the file was opened. */
  if (!fh)
    return 0;

  /* print the headings. */
  for (n = 0; n < 10; n++) fprintf(fh, " ");
  for (d = 0; d < D->nd; d++)
    fprintf(fh, "        Axis %2u", d + 1);

  /* print a newline. */
  fprintf(fh, "\n");

  /* print the points count. */
  fprintf(fh, "Points:   ");
  for (d = 0; d < D->nd; d++)
    fprintf(fh, "%15u", D->dims[d].sz);

  /* print a newline. */
  fprintf(fh, "\n");

  /* print the carrier frequencies. */
  fprintf(fh, "Obs (MHz):");
  for (d = 0; d < D->nd; d++)
    fprintf(fh, "%15.3f", D->dims[d].carrier);

  /* print a newline. */
  fprintf(fh, "\n");

  /* print the spectral widths. */
  fprintf(fh, "SW (Hz):  ");
  for (d = 0; d < D->nd; d++)
    fprintf(fh, "%15.3f", D->dims[d].width);

  /* print a newline. */
  fprintf(fh, "\n");

  /* print the spectral offsets. */
  fprintf(fh, "Off (Hz): ");
  for (d = 0; d < D->nd; d++)
    fprintf(fh, "%15.3f", D->dims[d].offset);

  /* print a newline. */
  fprintf(fh, "\n");

  /* print the real/complex modes. */
  fprintf(fh, "Mode:     ");
  for (d = 0; d < D->nd; d++)
    fprintf(fh, "%15s", D->dims[d].cx ? "Complex" : "Real");

  /* print a newline. */
  fprintf(fh, "\n");

  /* print the uniformity status. */
  fprintf(fh, "NUS:      ");
  for (d = 0; d < D->nd; d++)
    fprintf(fh, "%15s", D->dims[d].nus ? "True" : "False");

  /* print a newline. */
  fprintf(fh, "\n");

  /* print the nucleus name strings. */
  fprintf(fh, "Name:     ");
  for (d = 0; d < D->nd; d++)
    fprintf(fh, "%15s", D->dims[d].nuc);

  /* print a newline. */
  fprintf(fh, "\n");

  /* close the output file. */
  if (fname)
    fclose(fh);

  /* return success. */
  return 1;
}

/* datum_save(): saves the contents of an acquired NMR datum to a file.
 * @D: the datum to save data from.
 * @fname: the output filename.
 */
int datum_save (datum *D, const char *fname) {
  /* FIXME: implement datum_save() */

  /* return success. */
  return 1;
}

/* datum_load(): loads the contents of an acquired NMR datum from a file.
 * @D: the datum to load data from.
 * @fname: the input filename.
 */
int datum_load (datum *D, const char *fname) {
  /* FIXME: implement datum_load() */

  /* return success. */
  return 1;
}

/* datum_reorder_dims(): reorders the dimensions inside a datum structure
 * according to the initial ordering @order. the target ordering will be
 * achieved by sorting the values in @order until they are in increasing
 * order.
 *
 * @D: pointer to the target datum structure.
 * @order: initial dimension ordering.
 */
int datum_reorder_dims (datum *D, int *order) {
  /* declare required variables:
   * @i: insertion sort outer loop counter.
   * @j: insertion sort inner loop counter.
   * @dim: a swap location for dimension information.
   */
  unsigned int i, j;
  datum_dim dim;
  int *ord, swp;

  /* allocate a temporary array of indices. */
  ord = (int*) malloc(D->nd * sizeof(int));
  if (!ord)
    return 0;

  /* copy the ordering into the temporary array. */
  memcpy(ord, order, D->nd * sizeof(int));

  /* loop over the array of dimensions. */
  for (i = 1; i < D->nd; i++) {
    /* set the initial inner loop index. */
    j = i;

    /* loop over the unsorted dimensions. */
    while (j > 0 && ord[j - 1] > ord[j]) {
      /* swap the (j-1) and (j) elements. */
      dim = D->dims[j];
      D->dims[j] = D->dims[j - 1];
      D->dims[j - 1] = dim;

      /* swap the ordering of the indices. */
      swp = ord[j];
      ord[j] = ord[j - 1];
      ord[j - 1] = swp;

      /* decrement the inner loop counter. */
      j--;
    }
  }

  /* free the ordering array. */
  free(ord);

  /* return success. */
  return 1;
}

