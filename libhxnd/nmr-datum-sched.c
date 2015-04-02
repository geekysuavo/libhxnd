
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

/* datum_sched_read(): read a bruker/varian schedule file into the schedule
 * array of a datum structure.
 * @D: pointer to the datum structure to manipulate.
 * @fname: the input filename.
 */
int datum_sched_read (datum *D, const char *fname) {
  /* declare required variables:
   * @buf: character buffer for reading lines.
   * @tokv: array of token strings.
   * @tokc: number of strings in @tokc.
   * @d: number of schedule columns.
   * @n: number of schedule rows.
   * @i: general purpose loop index.
   * @sched: output schedule array.
   * @fh: input file handle.
   */
  char buf[N_BUF], **tokv;
  unsigned int toki, tokc;
  hx_index sched;
  int d, n, i;
  FILE *fh;

  /* initialize the results. */
  sched = NULL;
  d = n = 0;

  /* open the input file. */
  fh = fopen(fname, "rb");

  /* check that the input file was opened. */
  if (!fh)
    throw("failed to open '%s'", fname);

  /* initialize the array index. */
  i = 0;

  /* loop until we've read the entire file. */
  while (!feof(fh)) {
    /* read a new line from the file. */
    if (fgets(buf, N_BUF, fh)) {
      /* split the read line into tokens. */
      tokv = strsplit(buf, " ", &tokc);

      /* check that the split was successful. */
      if (!tokv || tokc < 1)
        throw("failed to split string '%s'", buf);

      /* trim the token strings. */
      strvtrim(tokv, tokc);

      /* store/check the number of array elements. */
      if (i == 0)
        d = tokc;
      else if (tokc != d)
        throw("unexpected token count %d", tokc);

      /* increment the row count. */
      n++;

      /* reallocate the schedule array. */
      sched = (hx_index) realloc(sched, d * n * sizeof(int));

      /* check that the reallocation succeeded. */
      if (!sched)
        throw("failed to reallocate schedule array");

      /* loop over the tokens. */
      for (toki = 0; toki < tokc; toki++, i++)
        sched[i] = atoi(tokv[toki]);

      /* free the string array. */
      strvfree(tokv, tokc);
    }
  }

  /* close the input file. */
  fclose(fh);

  /* store the identified array parameters. */
  D->d_sched = d;
  D->n_sched = n;

  /* store the constructed schedule array. */
  D->sched = sched;

  /* return success. */
  return 1;
}

/* datum_sched_free(): free the nonuniform schedule information contained
 * within a datum structure, if any exists.
 * @D: pointer to the datum to manipulate.
 */
int datum_sched_free (datum *D) {
  /* declare a required variable:
   * @d: datum dimension index.
   */
  int d;

  /* free the datum schedule array, if it exists. */
  if (D->sched)
    free(D->sched);

  /* re-initialize the datum schedule information. */
  D->sched = NULL;
  D->d_sched = 0;
  D->n_sched = 0;

  /* set all datum array dimensions to uniformly sampled. */
  for (d = 0; d < D->nd; d++)
    D->dims[d].nus = 0;

  /* return success. */
  return 1;
}

