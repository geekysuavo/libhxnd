
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
   * @idx: unpacked multidimensional array indices.
   * @pidx: packed linear array index.
   * @d: dimension loop counter.
   * @fh: output file handle.
   */
  unsigned int d;
  hx_index idx;
  int i, pidx;
  FILE *fh;
  real fi;

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
  idx = hx_index_alloc(D->array.k);

  /* check that allocation was successful. */
  if (!idx)
    throw("failed to allocate %d indices", D->array.k);

  /* iterate over the points in the core array. */
  do {
    /* pack the linear index. */
    hx_index_pack(D->array.k, D->array.sz, idx, &pidx);

    /* print the indices. */
    for (i = 0; i < D->array.k; i++) {
      fi = (real) idx[i] / (real) (D->dims[i].sz - 1) - 0.5;
      fprintf(fh, "%18.8e ", D->dims[i].offset + fi * D->dims[i].width);
    }

    /* print the coefficients. */
    for (i = 0; i < D->array.n; i++)
      fprintf(fh, "%18.8e ", D->array.x[i + D->array.n * pidx]);

    /* print a newline. */
    fprintf(fh, "\n");
  } while (hx_index_incr(D->array.k, D->array.sz, idx));

  /* close the output file. */
  fclose(fh);

  /* free the multidimensional index. */
  hx_index_free(idx);

  /* return success. */
  return 1;
}

