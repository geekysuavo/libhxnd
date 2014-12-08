
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

/* include the rnmrtk header. */
#include <hxnd/nmr-rnmrtk.h>

/* declare a buffer size for reading rnmrtk parameter files.
 */
#define N_BUF  256

/* rnmrtk_parfile(): build a parameter filename string.
 * @fname: the input data file name.
 */
char *rnmrtk_parfile (const char *fname) {
  /* declare a few required variables. */
  unsigned int n_fname;
  char *pfname;

  /* check that the string is non-null and has sufficient length. */
  if (!fname || strlen(fname) <= 4)
    return NULL;

  /* get the size of the new string. */
  n_fname = strlen(fname) + 2;

  /* allocate the parameter filename string. */
  pfname = (char*) malloc(n_fname * sizeof(char));

  /* check that allocation succeeded. */
  if (!pfname)
    return NULL;

  /* build the parameter filename string. */
  snprintf(pfname, n_fname, "%s", fname);
  pfname[strlen(pfname) - 4] = '.';
  pfname[strlen(pfname) - 3] = 'p';
  pfname[strlen(pfname) - 2] = 'a';
  pfname[strlen(pfname) - 1] = 'r';

  /* return the allocated filename. */
  return pfname;
}

/* rnmrtk_check_file(): check a file for the requisite properties that may
 * define it as an RNMRTK data file.
 * @fname: the input data file name.
 */
int rnmrtk_check_file (const char *fname) {
  /* declare a few required variables. */
  int have_dat, have_par;
  char *pfname;

  /* check if the data file exists. */
  have_dat = bytes_fexist(fname);

  /* build the parameter filename string. */
  pfname = rnmrtk_parfile(fname);

  /* check that the string was created successfully. */
  if (!pfname)
    return 0;

  /* check if the parameter file exists. */
  have_par = bytes_fexist(pfname);

  /* free the parameter filename string. */
  free(pfname);

  /* clear any errors that may have popped up. */
  traceback_clear();

  /* return whether both files exist. */
  return (have_dat && have_par);
}

/* rnmrtk_read_parms(): read the parameter file corresponding to an RNMRTK
 * data file.
 * @fname: the data filename.
 * @par: pointer to the output parameter structure.
 */
int rnmrtk_read_parms (const char *fname, struct rnmrtk_parms *par) {
  /* declare a few required variables:
   * @pfname: parameter file name.
   * @buf: buffer holding parameter file lines.
   * @nfields: number of strings in the fields array.
   * @fields: string array of space-delimited fields.
   * @i: general purpose loop counter.
   * @fh: parameter file handle.
   */
  char *pfname, buf[N_BUF];
  unsigned int i, nfields;
  char **fields;
  FILE *fh;

  /* declare variables required for parsing layout lines:
   * @dim: main dimension index.
   * @sub: sub-dimension index.
   * @pts: point count.
   */
  int dim, sub, pts;

  /* initialize the parameter structure contents. */
  memset(par, 0, sizeof(struct rnmrtk_parms));

  /* build the parameter filename string. */
  pfname = rnmrtk_parfile(fname);

  /* check that the string was created successfully. */
  if (!pfname)
    throw("failed to allocate parameter filename");

  /* open the parameter file. */
  fh = fopen(pfname, "rb");

  /* check that the file was opened. */
  if (!fh)
    throw("failed to open '%s'", pfname);

  /* loop until we've read the entire file. */
  while (!feof(fh)) {
    /* read a new line from the file. */
    if (fgets(buf, N_BUF, fh)) {
      /* trim trailing newlines from the string. */
      strnltrim((char*) buf);

      /* split the line by whitespace. */
      fields = strsplit(buf, " ", &nfields);

      /* trim the fields and convert everything to lowercase. */
      strvtrim(fields, nfields);
      strvcompact(fields, &nfields);
      strvtolower(fields, nfields);

      /* determine which parameter line we've tokenized. */
      if (strcmp(fields[0], RNMRTK_PARLINE_FORMAT) == 0) {
        /* parse all available format line fields. */
        if (nfields >= 2) {
          /* parse the endianness of the data. */
          if (strcmp(fields[1], RNMRTK_ENDIAN_BIG) == 0)
            par->endian = BYTES_ENDIAN_BIG;
          else if (strcmp(fields[1], RNMRTK_ENDIAN_LITTLE) == 0)
            par->endian = BYTES_ENDIAN_LITTLE;
          else
            throw("invalid endianness '%s'", fields[1]);
        }
        if (nfields >= 3) {
          /* parse the word style of the data. */
          if (strcmp(fields[2], RNMRTK_WTYPE_INT) == 0)
            par->isflt = 0;
          else if (strcmp(fields[2], RNMRTK_WTYPE_FLT) == 0)
            par->isflt = 1;
          else
            throw("invalid word type '%s'", fields[2]);
        }
        if (nfields >= 4) {
          /* parse the header size. */
          par->nheader = atoi(fields[3]);
        }
        if (nfields >= 5) {
          /* parse the record length. */
          par->reclen = atoi(fields[4]);
        }
        if (nfields >= 6) {
          /* parse the padding size. */
          par->nbegin = atoi(fields[5]);
        }
        if (nfields >= 7) {
          /* parse the end padding. */
          par->nend = atoi(fields[6]);
        }
      }
      else if (strcmp(fields[0], RNMRTK_PARLINE_DOM) == 0) {
        /* store the number of dimensions. */
        par->nd = nfields - 1;

        /* check the number of dimensions. */
        if (par->nd < 1 || par->nd > RNMRTK_MAXDIM)
          throw("invalid dimension count %u", par->nd);

        /* parse the 'dom' line fields. */
        for (i = 0; i < par->nd; i++) {
          /* check that the first character is 't'. */
          if (fields[i + 1][0] != 't' || strlen(fields[i + 1]) < 2)
            throw("invalid field '%s' on '%s'", fields[i + 1], fields[0]);

          /* parse the dimension index. */
          par->ord[i] = atoi(fields[i + 1] + 1);
        }
      }
      else if (strcmp(fields[0], RNMRTK_PARLINE_N) == 0) {
        /* check for correct field count. */
        if (nfields != 2 * par->nd + 1)
          throw("invalid field count of %u on '%s'", nfields, fields[0]);

        /* parse the 'n' line fields. */
        for (i = 0; i < par->nd; i++) {
          /* parse the size value. */
          par->sz[i] = atoi(fields[2 * i + 1]);

          /* parse the complex flag. */
          if (strcmp(fields[2 * i + 2], RNMRTK_NTYPE_REAL) == 0)
            par->cx[i] = 0;
          else if (strcmp(fields[2 * i + 2], RNMRTK_NTYPE_COMPLEX) == 0)
            par->cx[i] = 1;
          else
            throw("invalid real/complex field '%s'", fields[2 * i + 2]);
        }
      }
      else if (strcmp(fields[0], RNMRTK_PARLINE_LAYOUT) == 0) {
        /* parse the 'layout' line fields. */
        for (i = 1; i < nfields; i++) {
          /* attempt to parse each layout field form. */
          if (sscanf(fields[i], "t%d-%d:%d", &dim, &sub, &pts) == 3 &&
              dim > 0 && dim <= par->nd) {
            /* store the parsed point count. */
            par->layout[dim - 1][sub] = pts;
          }
          else if (sscanf(fields[i], "t%d:%d", &dim, &pts) == 2 &&
                   dim > 0 && dim <= par->nd) {
            /* store the parsed point count. */
            par->layout[dim - 1][1] = pts;
          }
          else
            throw("invalid %s field '%s'", fields[0], fields[i]);
        }
      }
      else if (strcmp(fields[0], RNMRTK_PARLINE_SF) == 0) {
        /* check for correct field count. */
        if (nfields != par->nd + 1)
          throw("invalid field count of %u on '%s'", nfields, fields[0]);

        /* parse the 'sf' line fields. */
        for (i = 0; i < par->nd; i++)
          par->sf[i] = atof(fields[i + 1]);
      }
      else if (strcmp(fields[0], RNMRTK_PARLINE_PPM) == 0) {
        /* check for correct field count. */
        if (nfields != par->nd + 1)
          throw("invalid field count of %u on '%s'", nfields, fields[0]);

        /* parse the 'ppm' line fields. */
        for (i = 0; i < par->nd; i++)
          par->ppm[i] = atof(fields[i + 1]);
      }
      else if (strcmp(fields[0], RNMRTK_PARLINE_QUAD) == 0) {
        /* check for correct field count. */
        if (nfields != par->nd + 1)
          throw("invalid field count of %u on '%s'", nfields, fields[0]);

        /* parse the 'quad' line fields. */
        for (i = 0; i < par->nd; i++) {
          /* determine the quadrature setting. */
          if (strcmp(fields[i + 1], RNMRTK_QUADSTR_TPPI) == 0)
            par->quad[i] = RNMRTK_QUAD_TPPI;
          else if (strcmp(fields[i + 1], RNMRTK_QUADSTR_STATES) == 0)
            par->quad[i] = RNMRTK_QUAD_STATES;
          else if (strcmp(fields[i + 1], RNMRTK_QUADSTR_STATESTPPI) == 0)
            par->quad[i] = RNMRTK_QUAD_STATESTPPI;
          else
            throw("invalid quadrature '%s'", fields[i + 1]);
        }
      }
      else if (strcmp(fields[0], RNMRTK_PARLINE_SW) == 0) {
        /* check for correct field count. */
        if (nfields != par->nd + 1)
          throw("invalid field count of %u on '%s'", nfields, fields[0]);

        /* parse the 'sw' line fields. */
        for (i = 0; i < par->nd; i++)
          par->sw[i] = atof(fields[i + 1]);
      }

      /* free the fields string array. */
      strvfree(fields, nfields);
    }
  }

  /* close the parameter file. */
  fclose(fh);

  /* free the parameter filename. */
  free(pfname);

  /* return success. */
  return 1;
}

/* rnmrtk_read(): read an RNMRTK data file into a real linear array.
 * @fname: the input data filename.
 * @x: the output array.
 */
int rnmrtk_read (const char *fname, hx_array *x) {
  /* declare a structure to hold parameter information. */
  struct rnmrtk_parms par;

  /* read the parameter file. */
  if (!rnmrtk_read_parms(fname, &par))
    throw("failed to read rnmrtk parameter file");

  /* FIXME: implement rnmrtk_read() */
  throw("unimplemented!");

  /* return success. */
  return 1;
}

/* rnmrtk_fill_datum(): load parameter information accompanying an RNMRTK
 * data file into an NMR datum structure.
 * @fname: the input data file name.
 * @D: pointer to the datum struct to fill.
 */
int rnmrtk_fill_datum (const char *fname, datum *D) {
  /* declare a structure to hold parameter information. */
  struct rnmrtk_parms par;

  /* read the parameter file. */
  if (!rnmrtk_read_parms(fname, &par))
    throw("failed to read rnmrtk parameter file");

  /* FIXME: implement rnmrtk_fill_datum() */
  throw("unimplemented!");

  /* return success. */
  return 1;
}

