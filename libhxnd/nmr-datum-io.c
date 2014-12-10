
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

/* declare a buffer size for reading bruker/varian schedule files.
 */
#define N_BUF  256

/* datum_load(): read a file (of any supported format) into a datum structure.
 * @D: pointer to the destination datum structure.
 * @fname: the input filename.
 */
int datum_load (datum *D, const char *fname) {
  /* initialize the datum contents. */
  datum_init(D);

  /* check the input filename. */
  if (!fname)
    throw("invalid input filename");

  /* determine the type of the datum. */
  D->type = datum_type_guess(fname);

  /* check that the type is supported. */
  if (D->type == DATUM_TYPE_UNDEFINED)
    throw("failed to identify format of '%s'", fname);

  /* decode the file parameters. */
  if (!datum_type_decode(D, fname, D->type))
    throw("failed to decode %s-format parameters from '%s'",
          datum_type_name(D->type), fname);

  /* read the array data from the file. */
  if (!datum_read_array(D))
    throw("failed to read %s-format array data from '%s'",
          datum_type_name(D->type), fname);

  /* return success. */
  return 1;
}

/* datum_print(): print the metadata associated with an acquired NMR datum.
 * @D: the datum to print data from.
 * @fname: the output filename.
 */
int datum_print (datum *D, const char *fname) {
  /* declare a few required variables:
   * @timedate: local calendar time structure.
   * @tmstr: date/time string.
   * @d: dimension loop counter.
   * @fh: the file handle used for writing.
   */
  struct tm *timedate;
  char tmstr[N_BUF];
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

  /* print the date line. */
  timedate = gmtime(&D->epoch);
  strftime(tmstr, N_BUF, "%c", timedate);
  fprintf(fh, "Date: %s\n", tmstr);

  /* print the array line. */
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
    }

    /* print a newline. */
    fprintf(fh, ")\n\n");
  }
  else
    fprintf(fh, "Unallocated.\n\n");

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
  fprintf(fh, "Domain:   ");
  for (d = 0; d < D->nd; d++)
    fprintf(fh, "%15s", D->dims[d].ft ? "Frequency" : "Time");

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

