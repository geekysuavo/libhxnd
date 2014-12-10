
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

/* define all possible parameter file statement strings supported by rnmrtk.
 */
#define RNMRTK_PARLINE_FORMAT  "format"
#define RNMRTK_PARLINE_DOM     "dom"
#define RNMRTK_PARLINE_N       "n"
#define RNMRTK_PARLINE_LAYOUT  "layout"
#define RNMRTK_PARLINE_SF      "sf"
#define RNMRTK_PARLINE_PPM     "ppm"
#define RNMRTK_PARLINE_QUAD    "quad"
#define RNMRTK_PARLINE_SW      "sw"

/* define all possible endianness strings in rnmrtk 'format' lines.
 */
#define RNMRTK_ENDIAN_BIG     "big-endian"
#define RNMRTK_ENDIAN_LITTLE  "little-endian"

/* define all possible word-type strings in rnmrtk 'format' lines.
 */
#define RNMRTK_WTYPE_INT  "int-32"
#define RNMRTK_WTYPE_FLT  "ieee-float"

/* define all possible quadrature strings in rnmrtk 'quad' lines.
 */
#define RNMRTK_QUADSTR_TPPI        "tppi"
#define RNMRTK_QUADSTR_STATES      "states"
#define RNMRTK_QUADSTR_STATESTPPI  "states-tppi"

/* define all possible real/complex strings in rnmrtk 'n' lines.
 */
#define RNMRTK_NTYPE_REAL     "r"
#define RNMRTK_NTYPE_COMPLEX  "c"

/* define the maximum number of dimensions supported by rnmrtk files.
 */
#define RNMRTK_MAXDIM   4
#define RNMRTK_MAXSUB  10

/* rnmrtk_quad: enumerated type describing the quadrature mode of a given
 * dimension in an RNMRTK data file.
 */
enum rnmrtk_quad {
  RNMRTK_QUAD_REAL       = 0x00,
  RNMRTK_QUAD_TPPI       = 0x01,
  RNMRTK_QUAD_STATES     = 0x02,
  RNMRTK_QUAD_STATESTPPI = 0x03
};

/* rnmrtk_parms: structure containing parsed parameters that correspond to
 * an RNMRTK data file.
 */
struct rnmrtk_parms {
  /* arguments parsed from parfile 'format' lines:
   * @endian: raw data byte ordering, little or big.
   * @isflt: whether the data is floating-point or integer.
   * @nheader: number of bytes in the data file header.
   * @reclen: number of real points per record.
   * @nbegin: number of bytes padding each data record.
   * @nend: number of bytes padding the end of the file.
   */
  enum byteorder endian;
  unsigned int isflt;
  unsigned int nheader;
  unsigned int reclen;
  unsigned int nbegin;
  unsigned int nend;

  /* arguments parsed from parfile 'dom' lines:
   * @ord: dimension ordering array, one-based.
   * @nd: dimension count.
   */
  int ord[RNMRTK_MAXDIM];
  unsigned int nd;

  /* arguments parsed from parfile 'n' lines:
   * @sz: array of sizes for each dimension.
   * @cx: array of complex flags for each dimension.
   */
  int sz[RNMRTK_MAXDIM], cx[RNMRTK_MAXDIM];

  /* arguments parsed from parfile 'layout' lines:
   * @layout: matrix of dimension,subdimension sizes.
   */
  int layout[RNMRTK_MAXDIM][RNMRTK_MAXSUB];

  /* arguments parsed from optional parfile lines:
   * @sf: spectrometer carrier frequencies.
   * @ppm: carrier offset values, in ppm.
   * @sw: spectral widths.
   */
  float sf[RNMRTK_MAXDIM];
  float ppm[RNMRTK_MAXDIM];
  float sw[RNMRTK_MAXDIM];
  enum rnmrtk_quad quad[RNMRTK_MAXDIM];
};

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

/* rnmrtk_guess(): check a file contains rnmrtk-format data.
 * @fname: the input data file name.
 */
int rnmrtk_guess (const char *fname) {
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

  /* return whether both files exist. */
  return (have_dat && have_par);
}

/* rnmrtk_decode(): read parameter information accompanying an rnmrtk
 * data file into a datum structure.
 * @D: pointer to the destination datum structure.
 * @fname: the input data file name.
 */
int rnmrtk_decode (datum *D, const char *fname) {
  /* @par: structure to hold parameter information.
   * @ord: dimension ordering array.
   * @d: dimension loop counter.
   */
  struct rnmrtk_parms par;
  int ord[RNMRTK_MAXDIM];
  unsigned int d;

  /* read the parameter file. */
  if (!rnmrtk_read_parms(fname, &par))
    throw("failed to read rnmrtk parameter file");

  /* store the dimensionality. */
  D->nd = par.nd;

  /* check the dimensionality. */
  if (D->nd < 1 || D->nd > RNMRTK_MAXDIM)
    throw("invalid dimensionality %u", D->nd);

  /* allocate the dimension parameter array. */
  D->dims = (datum_dim*) calloc(D->nd, sizeof(datum_dim));

  /* check that the dimension parameter array was allocated. */
  if (D->dims == NULL)
    throw("failed to allocate %u datum dimensions", D->nd);

  /* store the dimension information. */
  for (d = 0; d < D->nd; d++) {
    /* store the dimension size. */
    D->dims[d].td = D->dims[d].tdunif = par.layout[d][1];
    D->dims[d].sz = par.sz[d];

    /* store the complex flag. */
    if (par.cx[d])
      D->dims[d].cx = 1;

    /* store the spectral parameters. */
    D->dims[d].carrier = par.sf[d];
    D->dims[d].width = par.sw[d];
    D->dims[d].offset = par.ppm[d] / par.sf[d];

    /* determine quadrature information. */
    switch (par.quad[d]) {
      /* real, states: no extra processing. */
      case RNMRTK_QUAD_REAL:
      case RNMRTK_QUAD_STATES:
        break;

      /* tppi, states-tppi: sign alternation. */
      case RNMRTK_QUAD_TPPI:
      case RNMRTK_QUAD_STATESTPPI:
        D->dims[d].alt = 1;
        break;

      /* other: assume no extra processing. */
      default:
        break;
    }
  }

  /* build a local copy of the dimension ordering array. */
  for (d = 0; d < RNMRTK_MAXDIM; d++)
    ord[d] = (d < par.nd ? par.nd - par.ord[d] : par.nd);

  /* swap the ordering of the dimension parameters. */
  if (!datum_reorder_dims(D, ord))
    throw("failed to reorder dimensions");

  /* store the filename string. */
  D->fname = (char*) malloc((strlen(fname) + 1) * sizeof(char));
  if (D->fname)
    strcpy(D->fname, fname);

  /* store the datum type. */
  D->type = DATUM_TYPE_RNMRTK;
  D->endian = par.endian;

  /* return success. */
  return 1;
}

/* rnmrtk_encode(): write a datum structure in rnmrtk-format to a file.
 * @D: pointer to the source datum structure.
 * @fname: the output filename.
 */
int rnmrtk_encode (datum *D, const char *fname) {
  /* declare variables required for parameter output:
   * @par: structure to hold all data file parameters.
   */
  struct rnmrtk_parms par;

  /* FIXME: implement rnmrtk_encode() */

  /* return success. */
  return 1;
}

/* rnmrtk_array(): read an rnmrtk data file into a datum array.
 * @D: pointer to the source datum structure.
 */
int rnmrtk_array (datum *D) {
  /* declare a structure to hold parameter information. */
  struct rnmrtk_parms par;

  /* declare variables required for array raw input:
   * @wordsz: number of bytes per data word.
   * @isflt: whether the words are floats.
   * @ntrue: number of bytes in the input file.
   * @ncalc: estimated value of @ntrue.
   * @offhead: file header byte offset.
   * @offblk: block header byte offset.
   * @offend: block footer byte offset.
   * @nblks: number of data blocks.
   * @nwords: number of words per block.
   * @nalign: block alignment byte count.
   * @fh: input file handle.
   */
  unsigned int wordsz;
  unsigned int isflt;
  unsigned int ntrue;
  unsigned int ncalc;
  unsigned int offhead;
  unsigned int offblk;
  unsigned int offend;
  unsigned int nblks;
  unsigned int nwords;
  unsigned int nalign;
  FILE *fh;

  /* check that the input filename is valid. */
  if (D->fname == NULL)
    throw("invalid input filename");

  /* read the parameter file. */
  if (!rnmrtk_read_parms(D->fname, &par))
    throw("failed to read rnmrtk parfile for '%s'", D->fname);

  /* check the word type. */
  if (par.isflt) {
    /* set the float flag and the word size. */
    wordsz = sizeof(float);
    isflt = 1;
  }
  else {
    /* set the int flag and the word size. */
    wordsz = sizeof(int32_t);
    isflt = 0;
  }

  /* get the header and block byte offsets. */
  offhead = par.nheader;
  offblk = par.nbegin;
  offend = par.nend;

  /* get the total file size. */
  ntrue = bytes_size(D->fname);

  /* check that the file size is valid. */
  if (ntrue <= offhead)
    throw("invalid data file size of %u bytes", ntrue);

  /* get (or compute) the number of words per block. */
  if (par.reclen) {
    /* get the number from the parm structure. */
    nwords = par.reclen;

    /* determine how many blocks are in the file. */
    nblks = (ntrue - offhead) / (offblk + nwords * wordsz + offend);
  }
  else {
    /* assume a single data block. */
    nblks = 1;

    /* determine the size of the data block. */
    nwords = (ntrue - offhead - offblk - offend) / wordsz;
  }

  /* check that the computed blocking information matches the file size. */
  ncalc = offhead + nblks * (offblk + nwords * wordsz + offend);
  if (ncalc != ntrue)
    throw("expected file size %uB does not match actual %uB", ncalc, ntrue);

  /* check if block alignment is required. */
  if (offend) {
    /* build the block alignment size. */
    nalign = offblk + nwords * wordsz + offend;
  }
  else
    nalign = 0;

  /* open the input file for reading. */
  fh = fopen(D->fname, "rb");

  /* check that the file was opened. */
  if (!fh)
    throw("failed to open '%s'", D->fname);

  /* read data from the file into the output array. */
  if (!hx_array_fread_raw(fh, &D->array, par.endian, wordsz, isflt, offhead,
                          offblk, nblks, nwords, nalign))
    throw("failed to read raw data from '%s'", D->fname);

  /* close the input file. */
  fclose(fh);

  /* return success. */
  return 1;
}

