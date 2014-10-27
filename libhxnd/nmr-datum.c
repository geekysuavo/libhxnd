
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
 * @D: ponter to the datum to initialize.
 */
void datum_init (datum *D) {
  /* initialize the input data filename. */
  D->fname = NULL;

  /* initialize the datum type. */
  D->type = DATUM_TYPE_UNDEFINED;

  /* initialize the dimension count and dimensions array. */
  D->dims = NULL;
  D->nd = 0;

  /* indicate that the array is not allocated. */
  D->array_alloc = 0;
}

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
    throw("failed to open '%s'", fname);

  /* print the filename line. */
  fprintf (fh, "File:  %s\n", D->fname ? D->fname : "Unknown");

  /* print the first line. */
  fprintf(fh, "Array: ");
  if (D->array_alloc) {
    /* print the initial information. */
    fprintf(fh, "d = %d, k = %d, sz = (", D->array.d, D->array.k);

    /* loop over the dimensions. */
    for (d = 0; d < D->array.k; d++) {
      /* print the size. */
      fprintf(fh, "%d", D->array.sz[d]);

      /* print a delimiter if required. */
      if (d < D->array.k - 1)
        fprintf(fh, ", ");

      /* print an ending if required. */
      if (d == D->array.k - 1)
        fprintf(fh, ")\n");
    }

    /* print a newline. */
    fprintf(fh, "\n");
  }
  else
    fprintf(fh, "not allocated.\n");

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
    throw("failed to allocate %u indices", D->nd);

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

/* datum_refactor_array(): repacks, infills and deinterlaces the core array
 * structure of an NMR datum until it's dimensionality and complexity agree
 * with the dimension parameter values.
 * @D: pointer to the datum to manipulate.
 */
int datum_refactor_array (datum *D) {
  /* declare a few required variables:
   * @d: dimension loop counter.
   */
  unsigned int d;

  /* check that the array has been allocated. */
  if (!D->array_alloc)
    throw("array is unallocated");

  /* loop over the acquisition dimensions to refactor the nD array. */
  for (d = 0; d < D->nd; d++) {
    /* repack indirect dimensions in the array. */
    if (d > 0 && !hx_array_repack(&D->array, D->dims[d - 1].sz))
      throw("failed to repack array dimension %d", d);

    /* check if the current dimension is complex. */
    if (D->dims[d].cx) {
      /* de-interlace this dimension. */
      if (!hx_array_deinterlace(&D->array))
        throw("failed to deinterlace complex dimension %d", d);
    }
    else {
      /* increment the dimensionality without de-interlacing. */
      if (!hx_array_resize(&D->array,
            D->array.d + 1, D->array.k, D->array.sz))
        throw("failed to complex-promote dimension %d", d);
    }
  }

  /* return success. */
  return 1;
}

/* datum_read_array(): reads and refactors the array data of a datum structure
 * that has been initialized with a correct set of parameters.
 * @D: pointer to the datum to manipulate.
 */
int datum_read_array (datum *D) {
  /* declare a few required variables. */
  unsigned int d, nblk, szblk;

  /* check if the array has been allocated. */
  if (D->array_alloc)
    throw("array is already allocated");

  /* check that the filename is non-null and non-empty. */
  if (D->fname == NULL || strlen(D->fname) == 0)
    throw("filename is invalid");

  /* load based on the type of data. */
  if (D->type == DATUM_TYPE_BRUKER) {
    /* determine the data block size. */
    szblk = 4 * D->dims[0].td;

    /* determine the data block count. */
    for (d = 1, nblk = 1; d < D->nd; d++)
      nblk *= D->dims[d].td;

    /* check if the blocks are 1.0 KiB-aligned. */
    if (szblk % 1024 == 0) {
      /* yes. use a (faster) single read, because no gaps exist. */
      szblk *= nblk;
      nblk = 1;
    }

    /* load the raw data from the fid/ser file. */
    if (!bruker_read(D->fname, D->endian, nblk, szblk, &D->array))
      throw("failed to read bruker data");
  }
  else if (D->type == DATUM_TYPE_VARIAN) {
    /* load the raw data from the fid file. */
    if (!varian_read(D->fname, &D->array))
      throw("failed to read varian data");
  }
  else
    throw("unsupported data type %d", D->type);

  /* indicate that the array has been allocated. */
  D->array_alloc = 1;

  /* refactor the core array. */
  if (!datum_refactor_array(D))
    throw("failed to refactor array");

  /* return success. */
  return 1;
}

/* datum_free_array(): frees an allocated array structure from an NMR datum.
 * @D: pointer to the datum to manipulate.
 */
int datum_free_array (datum *D) {
  /* check if the array is allocated. */
  if (!D->array_alloc)
    return 1;

  /* free the data array. */
  hx_array_free(&D->array);

  /* indicate that the array has been de-allocated. */
  D->array_alloc = 0;

  /* return success. */
  return 1;
}

