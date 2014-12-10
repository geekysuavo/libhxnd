
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

/* include the text header. */
#include <hxnd/nmr-text.h>

/* text_encode(): write a datum structure in text-format to a file.
 * @D: pointer to the source structure.
 * @fname: the output filename.
 */
int text_encode (datum *D, const char *fname) {
  /* declare required variables:
   * @arr: unpacked multidimensional array indices.
   * @idx: packed linear array index.
   * @d: dimension loop counter.
   * @fh: output file handle.
   */
  int i, *arr, idx;
  unsigned int d;
  FILE *fh;

  /* open the output file. */
  if (fname)
    fh = fopen(fname, "wb");
  else
    fh = stdout;

  /* check that the file was opened. */
  if (!fh)
    throw("failed to open '%s'", fname);

  /* print heading information for each dimension. */
  for (d = 0; d < D->nd; d++) {
    /* print the dimension index and name. */
    fprintf(fh, "# Axis %2u ('%s'):\n", d + 1, D->dims[d].nuc);

    /* print the dimension sizes. */
    fprintf(fh, "# Points:   %15u\n", D->dims[d].sz);
    fprintf(fh, "# Total:    %15u\n", D->dims[d].td);

    /* print the spectral parameters. */
    fprintf(fh, "# Obs (MHz):%15.3f\n", D->dims[d].carrier);
    fprintf(fh, "# SW (Hz):  %15.3f\n", D->dims[d].width);
    fprintf(fh, "# Off (Hz): %15.3f\n", D->dims[d].offset);

    /* print a spacer. */
    fprintf(fh, "#\n");
  }

  /* allocate an index array for iteration. */
  arr = hx_array_index_alloc(D->array.k);

  /* check that allocation was successful. */
  if (!arr)
    throw("failed to allocate %d indices", D->array.k);

  /* iterate over the points in the core array. */
  idx = 0;
  do {
    /* print the indices. */
    for (i = 0; i < D->array.k; i++)
      fprintf(fh, "%6d ", arr[i]);

    /* print the coefficients. */
    for (i = 0; i < D->array.n; i++)
      fprintf(fh, "%18.8e ", D->array.x[i + D->array.n * idx]);

    /* print a newline. */
    fprintf(fh, "\n");

    /* increment the linear index. */
    idx++;
  } while (hx_array_index_incr(D->array.k, D->array.sz, arr));

  /* close the output file. */
  fclose(fh);

  /* free the index array. */
  free(arr);

  /* return success. */
  return 1;
}

