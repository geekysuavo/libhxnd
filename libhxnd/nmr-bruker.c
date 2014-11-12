
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

/* include the bruker header. */
#include <hxnd/nmr-bruker.h>

/* declare a buffer size for reading bruker parameter files.
 */
#define N_BUF  256

/* bruker_check_dir(): checks a directory for the requisite files for
 * loading bruker-formatted raw data.
 * @dname: the input directory name.
 */
int bruker_check_dir (const char *dname) {
  /* declare a few required variables. */
  int have_acqus, have_fid, have_ser;
  unsigned int n_fname;
  char *fname;
  FILE *fh;

  /* allocate a string for checking file existence. */
  n_fname = strlen(dname) + 16;
  fname = (char*) malloc(n_fname * sizeof(char));

  /* check that the string was allocated. */
  if (!fname)
    throw("failed to allocate %d-char buffer", n_fname);

  /* try to open the acqus file. */
  have_acqus = 0;
  snprintf(fname, n_fname, "%s/acqus", dname);
  fh = fopen(fname, "rb");

  /* check if the acqus file was opened. */
  if (fh) {
    /* yep. it exists. */
    have_acqus = 1;
    fclose(fh);
  }

  /* try to open the fid file. */
  have_fid = 0;
  snprintf(fname, n_fname, "%s/fid", dname);
  fh = fopen(fname, "rb");

  /* check if the fid file was opened. */
  if (fh) {
    /* yep. it exists. */
    have_fid = 1;
    fclose(fh);
  }

  /* try to open the ser file. */
  have_ser = 0;
  snprintf(fname, n_fname, "%s/ser", dname);
  fh = fopen(fname, "rb");

  /* check if the ser file was opened. */
  if (fh) {
    /* yep. it exists. */
    have_ser = 1;
    fclose(fh);
  }

  /* free the filename string. */
  free(fname);

  /* return success. */
  return (have_acqus && (have_fid || have_ser));
}

/* bruker_read_parms(): read any number of parameters from a bruker
 * 'acqu', 'acqus', 'proc' or 'procs' file. each parameter requested
 * must be provided as a type char, a key string, and a result pointer.
 * @fname: the parameter filename.
 * @n: the number of parameters to read.
 * @...: the parameters requested.
 */
int bruker_read_parms (const char *fname, unsigned int n, ...) {
  /* declare required variables for file/string parsing:
   * @buf: buffer string for each line of the file.
   * @attr: string array of parameter key/value pairs.
   * @ncpy: number of chars to copy for key/value pairs.
   * @nattr: length of @attr.
   * @fh: the input file handle.
   */
  char buf[N_BUF], **attr;
  unsigned int ncpy, nattr;
  FILE *fh;

  /* declare required variables for variable arguments parsing:
   * @i: parameter index, out of @n.
   * @nid: number of parameters parsed.
   * @vl: variable argument list structure.
   */
  unsigned int i, nid;
  va_list vl;

  /* declare required variables for key/value pair parsing:
   * @keys: key string array.
   * @typs: value type char array.
   * @vals: result pointer array.
   * @key: currently parsed key.
   * @val: currently parsed pointer.
   */
  char **keys, *typs;
  void **vals;
  char *key;
  void *val;

  /* open the file for reading. */
  fh = fopen(fname, "rb");

  /* check that the file was opened. */
  if (!fh)
    throw("failed to open '%s'", fname);

  /* initialize the variable arguments list. */
  va_start(vl, n);

  /* allocate memory to store the search options for each parameter. */
  keys = (char**) malloc(n * sizeof(char*));
  vals = (void**) malloc(n * sizeof(void*));
  typs = (char*) malloc(n * sizeof(char));

  /* check that the memory was allocated successfully. */
  if (!keys || !vals || !typs)
    throw("failed to allocate option arrays");

  /* loop through the expected number of parameters. */
  for (i = 0; i < n; i++) {
    /* read a type, a key and a pointer from the arguments list. */
    typs[i] = (char) va_arg(vl, int);
    key = va_arg(vl, char*);
    val = va_arg(vl, void*);

    /* allocate memory to hold the key string. */
    keys[i] = (char*) malloc((strlen(key) + 1) * sizeof(char));

    /* check that key string memory was allocated. */
    if (!keys[i])
      throw("failed to allocate key string %d", i);

    /* copy the key string and pointer into their arrays. */
    strcpy(keys[i], key);
    vals[i] = val;
  }

  /* free the variable arguments list and initialize the number
   * of read parameters.
   */
  va_end(vl);
  nid = 0;

  /* loop until we've read the entire file. */
  while (!feof(fh)) {
    /* read a new line from the file. */
    if (fgets(buf, N_BUF, fh)) {
      /* trim trailing newlines from the string. */
      strnltrim((char*) buf);

      /* check if the current line can contain a key/value pair. */
      if (strlen(buf) <= 3 || strncmp(buf, "##$", 3) || !strstr(buf, "= "))
        continue;

      /* split the string by an equals sign. */
      attr = strsplit(buf + 3, "= ", &nattr);

      /* loop through the parameters to search for. */
      for (i = 0; i < n; i++) {
        /* check if the current parameter key is a match. */
        if (strcmp(keys[i], attr[0]) == 0) {
          /* yep! increment the number of identified parameters. */
          nid++;

          /* act based on the parameter type. */
          switch (typs[i]) {
            /* integer. */
            case BRUKER_PARMTYPE_INT:
              *((int*) vals[i]) = atol(attr[1]);
              break;

            /* float */
            case BRUKER_PARMTYPE_FLOAT:
              *((float*) vals[i]) = atof(attr[1]);
              break;

            /* string */
            case BRUKER_PARMTYPE_STRING:
              ncpy = strlen(attr[1]) - 2;
              strncpy(vals[i], attr[1] + 1, ncpy);
              ((char*) vals[i])[ncpy] = '\0';
              break;

            /* other. */
            default:
              break;
          }
        }
      }

      /* free the string array. */
      strvfree(attr, nattr);
    }
  }

  /* close the input file. */
  fclose(fh);

  /* free each key string in the array. */
  for (i = 0; i < n; i++)
    free(keys[i]);

  /* free the temporary search option arrays. */
  free(keys);
  free(vals);
  free(typs);

  /* return the number of read parameters. */
  return nid;
}

/* bruker_read(): reads a bruker data file into a real linear array.
 * @fname: the input data filename.
 * @endianness: the data byte ordering.
 * @nblk: the number of blocks/fids.
 * @szblk: the size of block/fid.
 * @x: the output array.
 */
int bruker_read (const char *fname, enum byteorder endianness,
                 unsigned int nblk, unsigned int szblk,
                 hx_array *x) {
  /* declare a few required variables. */
  unsigned int n;
  uint8_t *bytes;

  /* read the data bytes in from the file. */
  bytes = bytes_read_bruker(fname, nblk, szblk, &n);

  /* check that the bytes were read successfully. */
  if (!bytes)
    throw("failed to read %u %u-byte blocks from '%s'", nblk, szblk, fname);

  /* build a real linear array from the byte data. */
  if (!bytes_toarray(bytes, n, endianness, 4, 0, x))
    throw("failed to convert bytes to array");

  /* free the read byte data. */
  free(bytes);

  /* return success. */
  return 1;
}

/* bruker_fix_grpdelay(): corrects a hypercomplex array loaded from a bruker
 * raw data file for group delay.
 * @x: pointer to the array to correct.
 * @grpdelay: group delay value, or -1 to use autodetection.
 */
int bruker_fix_grpdelay (hx_array *x, real grpdelay) {
  /* declare a few required variables:
   * @I: current absolute intensity value.
   * @Ire: current real intensity value.
   * @Iim: current imaginary intensity value.
   * @Imax: maximum absolute intensity value.
   * @gd: local variable used for group delay correction.
   * @gdmax: first trace index having maximal absolute intensity.
   */
  real I, Ire, Iim, Imax;
  int gd, gdmax;

  /* store the group delay value locally. */
  gd = (int) grpdelay;

  /* check if the group delay value needs autodetection. */
  if (gd == -1) {
    /* loop until the trace intensity maximum is located. */
    for (gd = 0, gdmax = 0, Imax = 0.0; gd < x->sz[0]; gd++) {
      /* get the real and imaginary components of the trace. */
      Ire = x->x[0 + gd * x->n];
      Iim = x->x[1 + gd * x->n];

      /* compute the trace magnitude. */
      I = sqrt(Ire * Ire + Iim * Iim);

      /* check if the current point has greater intensity than @Imax. */
      if (I > Imax) {
        /* store the new value. */
        gdmax = gd;
        Imax = I;
      }
    }

    /* store the autodetected group delay. */
    gd = gdmax;
  }

  /* left-shift each trace based on the group delay value. */
  if (!hx_array_shift(x, 0, -gd))
    throw("failed to correct bruker group delay");

  /* return success. */
  return 1;
}

/* bruker_fill_datum(): completely loads bruker raw data into an NMR datum
 * structure.
 * @dname: the input directory name.
 * @D: pointer to the datum struct to fill.
 */
int bruker_fill_datum (const char *dname, datum *D) {
  /* declare variables for filename generation:
   * @n_fname: buffer sizes of filename strings.
   * @fname_data: the 'fid' or 'ser' filename.
   * @fname_parm: the 'acqus' or 'acqu*s' filename.
   */
  unsigned int n_fname;
  char *fname_data;
  char *fname_parm;

  /* declare variables for identifying data dimensionality:
   * @acqus_parmode: the acqus 'PARMODE' parameter.
   * @d: dimension loop counter.
   * @ord: dimension index array for 'AQSEQ' corrections.
   */
  int acqus_parmode = -1;
  unsigned int d;
  int *ord;

  /* declare variables for parsing data byte ordering:
   * @endianness: the determined data byte ordering.
   * @acqus_bytorda: the acqus 'BYTORDA' parameter.
   */
  enum byteorder endianness;
  int acqus_bytorda = 0;

  /* declare variables for parsing dimension parameters:
   * @acqus_td: number of time-domain data points.
   * @acqus_nustd: number of 'nus' data points.
   * @acqus_aqmod: direct acquisition mode.
   * @acqus_fnmode: indirect acquisition mode.
   */
  real acqus_swh, acqus_sfo, acqus_offs;
  int acqus_td, acqus_nustd;
  enum bruker_aqseq acqus_aqseq;
  enum bruker_aqmod acqus_aqmod;
  enum bruker_fnmode acqus_fnmode;

  /* declare variables for correcting group delay:
   * @acqus_grpdly: group delay value.
   */
  real acqus_grpdly = -1.0;

  /* allocate memory for the fid/ser and acqu*s filenames. */
  n_fname = strlen(dname) + 16;
  fname_data = (char*) malloc(n_fname * sizeof(char));
  fname_parm = (char*) malloc(n_fname * sizeof(char));

  /* check that the filename strings were allocated. */
  if (!fname_data || !fname_parm)
    throw("failed to allocate filename strings");

  /* build the first parameter filename. */
  snprintf(fname_parm, n_fname, "%s/acqus", dname);

  /* parse an initial set of parameters from the acqus file. */
  if (bruker_read_parms(fname_parm, 2,
        BRUKER_PARMTYPE_INT, "PARMODE", &acqus_parmode,
        BRUKER_PARMTYPE_INT, "BYTORDA", &acqus_bytorda) != 2)
    throw("failed to get PARMODE/BYTORDA from '%s'", fname_parm);

  /* parse the group delay from the acqus file. */
  bruker_read_parms(fname_parm, 1,
    BRUKER_PARMTYPE_FLOAT, "GRPDLY", &acqus_grpdly);

  /* store the parsed group delay value. */
  D->grpdelay = acqus_grpdly;

  /* parse the acquisition sequence parameter, which will only succeed for
   * three-dimensional or higher data.
   */
  acqus_aqseq = 0;
  bruker_read_parms(fname_parm, 1,
    BRUKER_PARMTYPE_INT, "AQSEQ", &acqus_aqseq);

  /* compute the dimensionality of the data. */
  D->nd = acqus_parmode + 1;

  /* check the dimensionality. */
  if (D->nd < 1)
    throw("invalid dimensionality %u", D->nd);

  /* determine the byte ordering of the data. */
  if (acqus_bytorda == 0)
    endianness = BYTES_ENDIAN_LITTLE;
  else
    endianness = BYTES_ENDIAN_BIG;

  /* build the data filename. */
  if (D->nd > 1)
    snprintf(fname_data, n_fname, "%s/ser", dname);
  else
    snprintf(fname_data, n_fname, "%s/fid", dname);

  /* allocate the dimension parameter array. */
  D->dims = (datum_dim*) calloc(D->nd, sizeof(datum_dim));

  /* check that the dimension parameter array was allocated. */
  if (D->dims == NULL)
    throw("failed to allocate %u datum dimensions", D->nd);

  /* loop over the acquisition dimensions. */
  for (d = 0; d < D->nd; d++) {
    /* check if the parameter filename needs updating. */
    if (d)
      snprintf(fname_parm, n_fname, "%s/acqu%us", dname, d + 1);

    /* initialize the results. */
    acqus_td = acqus_nustd = 0;
    acqus_aqmod = -1;
    acqus_fnmode = -1;

    /* parse the important values from the parameter file. */
    if (bruker_read_parms(fname_parm, 1,
          BRUKER_PARMTYPE_INT, "TD", &acqus_td) != 1)
      throw("failed to get TD from '%s'", fname_parm);

    /* read any other (quasi-optional) parameters from the file. */
    bruker_read_parms(fname_parm, 7,
      BRUKER_PARMTYPE_INT, "NusTD", &acqus_nustd,
      BRUKER_PARMTYPE_INT, "AQ_mod", &acqus_aqmod,
      BRUKER_PARMTYPE_INT, "FnMODE", &acqus_fnmode,
      BRUKER_PARMTYPE_FLOAT, "SW_h", &acqus_swh,
      BRUKER_PARMTYPE_FLOAT, "O1", &acqus_offs,
      BRUKER_PARMTYPE_FLOAT, "SFO1", &acqus_sfo,
      BRUKER_PARMTYPE_STRING, "NUC1", D->dims[d].nuc);

    /* store the read parameters in the datum structure. */
    D->dims[d].td = D->dims[d].sz = acqus_td;
    D->dims[d].tdunif = acqus_nustd;

    /* store the read spectral parameters in the datum structure. */
    D->dims[d].carrier = acqus_sfo;
    D->dims[d].width = acqus_swh;
    D->dims[d].offset = acqus_offs;

    /* determine if the dimension is nonuniformly subsampled. */
    if (D->dims[d].tdunif && D->dims[d].td < D->dims[d].tdunif)
      D->dims[d].nus = 1;

    /* determine if the dimension is complex or real. */
    if ((d == 0 && acqus_aqmod != BRUKER_AQMOD_QF) ||
        (d > 0 && acqus_fnmode != BRUKER_FNMODE_QF)) {
      /* set the complexity flag and correct the actual points count. */
      D->dims[d].cx = 1;
      D->dims[d].sz /= 2;
    }
  }

  /* check if the indirect dimensions were acquired in reverse order. */
  if (D->nd >= 3 && acqus_aqseq != BRUKER_AQSEQ_321) {
    /* allocate an array of dimension indices. */
    ord = (int*) calloc(D->nd, sizeof(int));
    if (!ord)
      throw("failed to allocate %u indices", D->nd);

    /* store the reversed dimension order implied by BRUKER_AQSEQ_312 */
    for (d = 1, ord[0] = 0; d < D->nd; d++)
      ord[d] = D->nd - d;

    /* swap the ordering of the dimension parameters. */
    if (!datum_reorder_dims(D, ord))
      throw("failed to reorder dimensions");

    /* free the dimension index array. */
    free(ord);
  }

  /* free the allocated filename strings. */
  D->fname = fname_data;
  free(fname_parm);

  /* store the datum type. */
  D->type = DATUM_TYPE_BRUKER;
  D->endian = endianness;

  /* return success. */
  return 1;
}

