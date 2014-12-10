
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

/* define constant parameter type characters for acqus file parsing.
 */
#define BRUKER_PARMTYPE_INT     'i'
#define BRUKER_PARMTYPE_FLOAT   'f'
#define BRUKER_PARMTYPE_STRING  's'

/* bruker_aqmod: enumerated type for direct acquisition modes.
 */
enum bruker_aqmod {
  BRUKER_AQMOD_QF,    /* single-channel (real) detection. */
  BRUKER_AQMOD_QSIM,  /* simultaneous quadrature detection. */
  BRUKER_AQMOD_QSEQ,  /* sequential quadrature detection. */
  BRUKER_AQMOD_DQD    /* digital quadrature detection. */
};

/* bruker_fnmode: enumerated type for indirect acquisition modes.
 */
enum bruker_fnmode {
  BRUKER_FNMODE_UNDEFINED,
  BRUKER_FNMODE_QF,
  BRUKER_FNMODE_QSEQ,
  BRUKER_FNMODE_TPPI,
  BRUKER_FNMODE_STATES,
  BRUKER_FNMODE_STATESTPPI,
  BRUKER_FNMODE_GRADIENT
};

/* bruker_aqseq: (3+)-dimensional acquisition ordering scheme.
 */
enum bruker_aqseq {
  BRUKER_AQSEQ_321,
  BRUKER_AQSEQ_312
};

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

/* bruker_guess(): check whether a directory contains bruker-format data.
 * @dname: the input directory name.
 */
int bruker_guess (const char *dname) {
  /* declare a few required variables. */
  int have_acqus, have_fid, have_ser;
  unsigned int n_fname;
  char *fname;

  /* allocate a string for checking file existence. */
  n_fname = strlen(dname) + 16;
  fname = (char*) malloc(n_fname * sizeof(char));

  /* check that the string was allocated. */
  if (!fname)
    throw("failed to allocate %u-char buffer", n_fname);

  /* check if the acqus file exists. */
  snprintf(fname, n_fname, "%s/acqus", dname);
  have_acqus = bytes_fexist(fname);

  /* check if the fid file exists. */
  snprintf(fname, n_fname, "%s/fid", dname);
  have_fid = bytes_fexist(fname);

  /* check if the ser file exists. */
  snprintf(fname, n_fname, "%s/ser", dname);
  have_ser = bytes_fexist(fname);

  /* free the filename string. */
  free(fname);

  /* return the result. */
  return (have_acqus && (have_fid || have_ser));
}

/* bruker_decode(): decode bruker parameter data into a datum structure.
 * @D: pointer to the destination datum structure.
 * @dname: the input directory name.
 */
int bruker_decode (datum *D, const char *dname) {
  /* declare variables for filename generation:
   * @n_fname: buffer sizes of filename strings.
   * @fname_data: the 'fid' or 'ser' filename.
   * @fname_parm: the 'acqus' or 'acqu*s' filename.
   */
  unsigned int n_fname;
  char *fname_data;
  char *fname_parm;

  /* declare a variable for date/time extraction:
   * @acqus_date: seconds since the epoch of data collection.
   */
  int acqus_date = 0;

  /* declare variables for identifying data dimensionality:
   * @acqus_parmode: the acqus 'PARMODE' parameter.
   * @d: dimension loop counter.
   * @ord: dimension index array for 'AQSEQ' corrections.
   */
  int acqus_parmode = -1;
  unsigned int d, nd;
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

  /* parse the date from the acqus file. */
  bruker_read_parms(fname_parm, 1,
    BRUKER_PARMTYPE_INT, "DATE", &acqus_date);

  /* store the parsed date value. */
  D->epoch = (time_t) acqus_date;

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
  nd = acqus_parmode + 1;

  /* check the dimensionality. */
  if (nd < 1)
    throw("invalid dimensionality %u", nd);

  /* determine the byte ordering of the data. */
  if (acqus_bytorda == 0)
    endianness = BYTES_ENDIAN_LITTLE;
  else
    endianness = BYTES_ENDIAN_BIG;

  /* build the data filename. */
  if (nd > 1)
    snprintf(fname_data, n_fname, "%s/ser", dname);
  else
    snprintf(fname_data, n_fname, "%s/fid", dname);

  /* allocate the dimension parameter array. */
  if (!datum_realloc_dims(D, nd))
    throw("failed to allocate dimension array");

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
    if (D->dims[d].tdunif && D->dims[d].td != D->dims[d].tdunif)
      D->dims[d].nus = 1;

    /* determine if the dimension is complex or real. */
    if ((d == 0 && acqus_aqmod != BRUKER_AQMOD_QF) ||
        (d > 0 && acqus_fnmode != BRUKER_FNMODE_QF)) {
      /* set the complexity flag and correct the actual points count. */
      D->dims[d].cx = 1;
      D->dims[d].sz /= 2;
    }

    /* determine if the dimension requires sign alternation. */
    if (d > 0 &&
        (acqus_fnmode == BRUKER_FNMODE_QSEQ ||
         acqus_fnmode == BRUKER_FNMODE_TPPI ||
         acqus_fnmode == BRUKER_FNMODE_STATESTPPI))
      D->dims[d].alt = 1;

    /* determine if the dimension requires gradient arithmetic. */
    if (d > 0 && acqus_fnmode == BRUKER_FNMODE_GRADIENT)
      D->dims[d].genh = 1;
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

  /* build the schedule filename. */
  snprintf(fname_parm, n_fname, "%s/nuslist", dname);

  /* read the schedule file, if one exists. */
  if (bytes_fexist(fname_parm) && !datum_read_sched(D, fname_parm))
    throw("failed to read bruker schedule");

  /* free the allocated filename strings. */
  D->fname = fname_data;
  free(fname_parm);

  /* store the datum type. */
  D->type = DATUM_TYPE_BRUKER;
  D->endian = endianness;

  /* return success. */
  return 1;
}

/* bruker_array(): read a bruker data file into a datum array.
 * @D: the destination datum structure.
 */
int bruker_array (datum *D) {
  /* declare a few required variables:
   * @d: dimension loop counter.
   * @szblk: the number of words per data block.
   * @nblk: the number of data blocks.
   * @fh: the input data file handle.
   */
  unsigned int d, szblk, nblk;
  FILE *fh;

  /* determine the block size. */
  szblk = D->dims[0].td;

  /* determine the block count. */
  for (d = 1, nblk = 1; d < D->nd; d++)
    nblk *= D->dims[d].td;

  /* check that the input filename is valid. */
  if (D->fname == NULL)
    throw("invalid input filename");

  /* open the input file for reading. */
  fh = fopen(D->fname, "rb");

  /* check that the file was opened. */
  if (!fh)
    throw("failed to open '%s'", D->fname);

  /* read data from the file into the output array. */
  if (!hx_array_fread_raw(fh, &D->array, D->endian, 4, 0,
                          0, 0, nblk, szblk, 1024))
    throw("failed to read raw data from '%s'", D->fname);

  /* close the input file. */
  fclose(fh);

  /* return success. */
  return 1;
}

/* bruker_post(): correct a hypercomplex array loaded from a bruker raw data
 * file for group delay.
 * @D: pointer to the datum structure.
 */
int bruker_post (datum *D) {
  /* declare a few required variables:
   * @I: current absolute intensity value.
   * @Ire: current real intensity value.
   * @Iim: current imaginary intensity value.
   * @Imax: maximum absolute intensity value.
   * @gd: local variable used for group delay correction.
   * @gdmax: first trace index having maximal absolute intensity.
   * @sznew: array of new size values.
   * @i: dimension loop counter.
   * @x: datum array structure pointer.
   */
  real I, Ire, Iim, Imax;
  int i, gd, gdmax, *sznew;
  hx_array *x;

  /* store the group delay value and array pointer locally. */
  gd = (int) D->grpdelay;
  x = &D->array;

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

  /* allocate the new size array. */
  sznew = hx_array_index_alloc(x->k);

  /* check that the size array was allocated. */
  if (!sznew)
    throw("failed to allocate %d indices", x->k);

  /* fill the new size array. */
  for (i = 0; i < x->k; i++)
    sznew[i] = x->sz[i];

  /* left-shift each trace based on the group delay value. */
  if (!hx_array_shift(x, 0, -gd))
    throw("failed to correct bruker group delay");

  /* correct the first dimension of the size array. */
  sznew[0] -= gd;

  /* crop the shifted array. */
  if (!hx_array_resize(x, x->d, x->k, sznew))
    throw("failed to crop bruker group delay points");

  /* free the allocated size array. */
  free(sznew);

  /* correct the datum fields. */
  D->dims[0].sz -= gd;
  D->grpdelay = 0.0;

  /* return success. */
  return 1;
}

