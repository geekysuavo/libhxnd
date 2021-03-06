
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

/* include the varian header. */
#include <hxnd/nmr-varian.h>

/* declare a buffer size for reading varian parameter files.
 */
#define N_BUF  256

/* define constant parameter type characters for procpar file parsing.
 */
#define VARIAN_PARMTYPE_INT     'i'
#define VARIAN_PARMTYPE_INTS    'I'
#define VARIAN_PARMTYPE_FLOAT   'f'
#define VARIAN_PARMTYPE_STRING  's'

/* define bit masks for file header status words.
 */
#define VARIAN_HDR_S_DATA     0x0001
#define VARIAN_HDR_S_SPEC     0x0002
#define VARIAN_HDR_S_32       0x0004
#define VARIAN_HDR_S_FLOAT    0x0008
#define VARIAN_HDR_S_COMPLEX  0x0010
#define VARIAN_HDR_S_HYPERCX  0x0020
#define VARIAN_HDR_S_ACQPAR   0x0080
#define VARIAN_HDR_S_SECND    0x0100
#define VARIAN_HDR_S_TRANSF   0x0200
#define VARIAN_HDR_S_NP       0x0800
#define VARIAN_HDR_S_NF       0x1000
#define VARIAN_HDR_S_NI       0x2000
#define VARIAN_HDR_S_NI2      0x4000

/* varian_hdr_file: structure definition for varian file headers.
 */
struct varian_hdr_file {
  uint32_t nblocks;  /* number of data blocks. */
  uint32_t ntraces;  /* number of traces per block. */
  uint32_t np;       /* number of elements per trace. */
  uint32_t ebytes;   /* number of bytes per element. */
  uint32_t tbytes;   /* number of bytes per trace. */
  uint32_t bbytes;   /* number of bytes per block. */
  uint16_t vers_id;  /* software version, file_id status bits. */
  uint16_t status;   /* status of whole file. */
  uint32_t nheaders; /* number of block headers per block. */
};

/* varian_hdr_blk: structure definition for varian block headers.
 */
struct varian_hdr_blk {
  uint16_t scale;    /* scaling factor. */
  uint16_t status;   /* status of data in block. */
  uint16_t index;    /* block index. */
  uint16_t mode;     /* mode of data in block. */
  uint32_t ctcount;  /* ct value for fid. */
  float lpval;       /* f2 (2D-f1) left phase in phasefile. */
  float rpval;       /* f2 (2D-f1) right phase in phasefile. */
  float lvl;         /* level drift compensation. */
  float tlt;         /* tilt drift compensation. */
};

/* varian_hdr_ext: structure definition for varian extended block headers.
 */
struct varian_hdr_ext {
  uint16_t s_spare1;
  uint16_t status;    /* status word for block header. */
  uint16_t s_spare2;
  uint16_t s_spare3;
  uint32_t l_spare1;
  float lpval1;       /* 2D-f2 left phase. */
  float rpval1;       /* 2D-f2 right phase. */
  float f_spare1;
  float f_spare2;
};

/* varian_read_parms(): read any number of parameters from a varian
 * 'procpar' file. each parameter requested must be provided as a
 * type char, a key string, and a result pointer.
 * @fname: the parameter filename.
 * @n: the number of parameters to read.
 * @...: the parameters requested.
 */
int varian_read_parms (const char *fname, unsigned int n, ...) {
  /* declare required variables for file/string parsing:
   * @buf: buffer string for each line of the file.
   * @attr: string array of parameter attributes.
   * @fields: string array of parameter values.
   * @ncpy: number of chars to copy for strings.
   * @nattr: length of @attr.
   * @nfields: length of @fields.
   * @fh: the input file handle.
   */
  char buf[N_BUF], **attr, **fields;
  unsigned int ncpy, nattr, nfields;
  FILE *fh;

  /* declare required variables for variable arguments parsing:
   * @i: parameter index, out of @n.
   * @nid: number of parameters parsed.
   * @vl: variable argument list structure.
   */
  unsigned int i, j, nid;
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

      /* split the line by whitespace. */
      attr = strsplit(buf, " ", &nattr);

      /* check if we've hit an attributes line. */
      if (nattr == 11) {
        /* loop through the parameters to search for. */
        for (i = 0; i < n; i++) {
          /* check if the current parameter key is a match. */
          if (strcmp(keys[i], attr[0]) == 0) {
            /* yep! increment the number of identified parameters. */
            nid++;

            /* tokenize the next line. */
            if (!fgets(buf, N_BUF, fh))
              continue;

            /* trim trailing newlines from the string. */
            strnltrim((char*) buf);

            /* split the next line by whitespace or quotes. */
            if (typs[i] == VARIAN_PARMTYPE_STRING)
              fields = strsplit(buf, "\"", &nfields);
            else
              fields = strsplit(buf, " ", &nfields);

            /* act based on the parameter type. */
            switch (typs[i]) {
              /* integer. */
              case VARIAN_PARMTYPE_INT:
                *((int*) vals[i]) = atol(fields[1]);
                break;

              /* integer array. */
              case VARIAN_PARMTYPE_INTS:
                for (j = 0; j < nfields && strlen(fields[j]); j++)
                  ((int*) vals[i])[j] = atol(fields[j]);
                break;

              /* float */
              case VARIAN_PARMTYPE_FLOAT:
                *((float*) vals[i]) = atof(fields[1]);
                break;

              /* string */
              case VARIAN_PARMTYPE_STRING:
                ncpy = strlen(fields[1]);
                strncpy(vals[i], fields[1], ncpy);
                ((char*) vals[i])[ncpy] = '\0';
                break;

              /* other. */
              default:
                break;
            }

            /* free the fields string array. */
            strvfree(fields, nfields);
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

/* varian_read_header(): reads a varian data file header into a struct.
 * @fname: the data filename.
 * @endianness: byte order result pointer.
 * @hdr: header result pointer.
 */
int varian_read_header (const char *fname,
                        enum byteorder *endianness,
                        struct varian_hdr_file *hdr) {
  /* declare a few required variables. */
  uint8_t *bytes;

  /* read in the file header bytes. */
  bytes = bytes_read_block(fname, 0, sizeof(struct varian_hdr_file));

  /* check that the bytes were read successfully. */
  if (!bytes)
    throw("failed to read file header from '%s'", fname);

  /* copy the header bytes onto the header structure. */
  memcpy(hdr, bytes, sizeof(struct varian_hdr_file));

  /* check if the bytes per element value is exceedingly large. if so,
   * it's highly probable that the read bytes were ordered in reverse of
   * the current machine usage and need to be swapped.
   */
  if (hdr->ebytes > 4096) {
    /* byte-swap the first six 32-bit words. */
    bytes_swap_u32(&hdr->nblocks);
    bytes_swap_u32(&hdr->ntraces);
    bytes_swap_u32(&hdr->np);
    bytes_swap_u32(&hdr->ebytes);
    bytes_swap_u32(&hdr->tbytes);
    bytes_swap_u32(&hdr->bbytes);

    /* byte-swap the next two 16-bit words. */
    bytes_swap_u16(&hdr->vers_id);
    bytes_swap_u16(&hdr->status);

    /* byte-swap the last 32-bit word. */
    bytes_swap_u32(&hdr->nheaders);

    /* opposite endianness. */
    *endianness = bytes_get_nonnative();
  }
  else {
    /* same endianness. */
    *endianness = bytes_get_native();
  }

  /* free the read bytes. */
  free(bytes);

  /* return success. */
  return 1;
}

/* varian_guess(): check whether a directory contains varian-format raw data.
 * @dname: the input directory name.
 */
int varian_guess (const char *dname) {
  /* declare a few required variables. */
  int have_procpar, have_fid;
  unsigned int n_fname;
  char *fname;

  /* allocate a string for checking file existence. */
  n_fname = strlen(dname) + 16;
  fname = (char*) malloc(n_fname * sizeof(char));

  /* check that the string was allocated. */
  if (!fname)
    throw("failed to allocate %u-char buffer", n_fname);

  /* check if the procpar file exists. */
  snprintf(fname, n_fname, "%s/procpar", dname);
  have_procpar = bytes_fexist(fname);

  /* check if the fid file exists.. */
  snprintf(fname, n_fname, "%s/fid", dname);
  have_fid = bytes_fexist(fname);

  /* free the filename string. */
  free(fname);

  /* return the result. */
  return (have_procpar && have_fid);
}

/* varian_array(): read a varian raw data file into a datum array.
 * @D: the datum structure to populate.
 */
int varian_array (datum *D) {
  /* declare a few required variables.
   */
  unsigned int isflt, nblk, szblk, offblk, offhead;
  struct varian_hdr_file hdr;
  enum byteorder endian;
  FILE *fh;

  /* check that the filename is valid. */
  if (D->fname == NULL)
    throw("invalid input filename");

  /* read the file header. */
  if (!varian_read_header(D->fname, &endian, &hdr))
    throw("failed to read header from '%s'", D->fname);

  /* compute block and header sizes for data loading. */
  nblk = hdr.nblocks;
  szblk = hdr.ntraces * hdr.tbytes;
  offblk = hdr.bbytes - szblk;
  offhead = sizeof(struct varian_hdr_file);

  /* determine whether the data is floating point or not. */
  isflt = (hdr.status & VARIAN_HDR_S_FLOAT ? 1 : 0);

  /* open the input file for reading. */
  fh = fopen(D->fname, "rb");

  /* check that the file was opened. */
  if (!fh)
    throw("failed to open '%s'", D->fname);

  /* read data from the file into the output array. */
  if (!hx_array_fread_raw(fh, &D->array, endian, hdr.ebytes, isflt,
                          offhead, offblk, nblk,
                          szblk / hdr.ebytes, 0))
    throw("failed to read raw data from '%s'", D->fname);

  /* close the input file. */
  fclose(fh);

  /* return success. */
  return 1;
}

/* varian_count_dims(): counts the number of nonzero-size dimensions in
 * a varian procpar file.
 * @fname: the input procpar filename.
 */
unsigned int varian_count_dims (const char *fname) {
  /* declare a few required variables:
   * @d: dimension loop counter.
   * @n: parameter count.
   */
  unsigned int d, n = 0;
  char parmstr[8];
  int v;

  /* try to find the first dimension point count. */
  v = 0;
  if (varian_read_parms(fname, 1,
        VARIAN_PARMTYPE_INT, "np", &v) == 1 &&
        v > 1)
    n++;
  else
    return 0;

  /* try to find the second dimension point count. */
  v = 0;
  if (varian_read_parms(fname, 1,
        VARIAN_PARMTYPE_INT, "ni", &v) == 1 &&
        v > 1)
    n++;
  else
    return n;

  /* loop over the range of supported dimensions. */
  for (d = 2; d < 32; d++) {
    /* build the parameter name string. */
    snprintf(parmstr, 8, "ni%u", d);

    /* try to find the point count. */
    v = 0;
    if (varian_read_parms(fname, 1,
          VARIAN_PARMTYPE_INT, parmstr, &v) == 1 &&
          v > 1)
      n++;
    else
      return n;
  }

  /* return the computed counter. */
  return n;
}

/* varian_parse_date(): computes the number of seconds past the epoch of
 * varian data collection from the date string in a procpar file.
 * @fname: the input procpar filename.
 */
time_t varian_parse_date (const char *fname) {
  /* declare a few required variables:
   */
  char tstr[N_BUF];
  struct tm ts;
  time_t t;

  /* initialize the time structure. */
  memset(&ts, 0, sizeof(struct tm));
  t = 0;

  /* parse the date string from the procpar file. */
  varian_read_parms(fname, 1, VARIAN_PARMTYPE_STRING, "date", tstr);

  /* attempt to parse the month from the date string. */
  ts.tm_mon = (
    strncmp(tstr, "Jan", 3) == 0 ?  0 :
    strncmp(tstr, "Feb", 3) == 0 ?  1 :
    strncmp(tstr, "Mar", 3) == 0 ?  2 :
    strncmp(tstr, "Apr", 3) == 0 ?  3 :
    strncmp(tstr, "May", 3) == 0 ?  4 :
    strncmp(tstr, "Jun", 3) == 0 ?  5 :
    strncmp(tstr, "Jul", 3) == 0 ?  6 :
    strncmp(tstr, "Aug", 3) == 0 ?  7 :
    strncmp(tstr, "Sep", 3) == 0 ?  8 :
    strncmp(tstr, "Oct", 3) == 0 ?  9 :
    strncmp(tstr, "Nov", 3) == 0 ? 10 :
    strncmp(tstr, "Dec", 3) == 0 ? 11 : 0);

  /* attempt to parse the numeric components from the date string. */
  sscanf(tstr, "%*s %d %d", &ts.tm_mday, &ts.tm_year);

  /* subtract 1900 from the year structure field. */
  ts.tm_year -= 1900;

  /* convert the time structure values into an epoch offset. */
  t = mktime(&ts);

  /* return the computed time. */
  return t;
}

/* varian_decode(): decode varian acquisition parameters into a datum
 * structure.
 * @D: pointer to the destination datum structure.
 * @dname: the input directory name.
 */
int varian_decode (datum *D, const char *dname) {
  /* declare variables for filename generation:
   * @n_fname: buffer sizes of filename strings.
   * @fname_data: the 'fid' filename string.
   * @fname_parm: the 'procpar' filename string.
   */
  unsigned int n_fname;
  char *fname_data;
  char *fname_parm;

  /* declare variables for looping over dimensions:
   * @d: dimension loop counter (integer).
   * @nd: dimension count.
   * @dstr: dimension loop counter (string).
   */
  unsigned int d, nd;
  char dstr[8];

  /* declare variables for parsing dimension parameters:
   * @parm_npmax: number of Nyquist points for NLS data.
   * @parm_np: number of sampled points.
   * @parmstr: parameter string.
   */
  real parm_swh, parm_sfrq, parm_rfp;
  int parm_np, parm_npmax, parm_phase[8];
  char parmstr[N_BUF];

  /* declare variables for parsing dimension ordering:
   * @arry: string array of 'array' parameter.
   * @narry: number of strings in @arry.
   */
  unsigned int i, narry;
  char **arry;
  int *ord;

  /* allocate memory for the fid and procpar filenames. */
  n_fname = strlen(dname) + 16;
  fname_data = (char*) malloc(n_fname * sizeof(char));
  fname_parm = (char*) malloc(n_fname * sizeof(char));

  /* check that the filename strings were allocated. */
  if (!fname_data || !fname_parm)
    throw("failed to allocate filename strings");

  /* build the procpar filename. */
  snprintf(fname_parm, n_fname, "%s/procpar", dname);

  /* build the fid filename. */
  snprintf(fname_data, n_fname, "%s/fid", dname);

  /* store the parsed date value. */
  D->epoch = varian_parse_date(fname_parm);;

  /* algorithmically guess the dimensionality of the data. */
  nd = varian_count_dims(fname_parm);

  /* check the dimensionality. */
  if (nd < 1)
    throw("invalid dimensionality %u", nd);

  /* allocate the dimension parameter array. */
  if (!datum_dims_realloc(D, nd))
    throw("failed to allocate dimension array");

  /* loop over the acquisition dimensions. */
  for (d = 0; d < D->nd; d++) {
    /* build a string version of the dimension index. */
    snprintf(dstr, 8, "%u", d);

    /* parse the number of points. */
    snprintf(parmstr, N_BUF, "n%c%s", d > 0 ? 'i' : 'p', d > 1 ? dstr : "");
    if (varian_read_parms(fname_parm, 1,
          VARIAN_PARMTYPE_INT, parmstr, &parm_np) != 1)
      throw("failed to get %s from '%s'", parmstr, fname_parm);

    /* parse the maximum number of points. */
    parm_npmax = 0;
    snprintf(parmstr, N_BUF, "n%cmax%s",
             d > 0 ? 'i' : 'p',
             d > 1 ? dstr : "");
    varian_read_parms(fname_parm, 1,
      VARIAN_PARMTYPE_INT, parmstr, &parm_npmax);

    /* parse the carrier base frequency. */
    snprintf(parmstr, N_BUF, "%cfrq%s", d > 0 ? 'd' : 's', d > 1 ? dstr : "");
    if (varian_read_parms(fname_parm, 1,
          VARIAN_PARMTYPE_FLOAT, parmstr, &parm_sfrq) != 1)
      throw("failed to get %s from '%s'", parmstr, fname_parm);

    /* parse the spectral width. */
    snprintf(parmstr, N_BUF, "sw%s", d > 0 ? dstr : "");
    if (varian_read_parms(fname_parm, 1,
          VARIAN_PARMTYPE_FLOAT, parmstr, &parm_swh) != 1)
      throw("failed to get %s from '%s'", parmstr, fname_parm);

    /* parse the spectral offset. */
    snprintf(parmstr, N_BUF, "rfp%s", d > 0 ? dstr : "");
    if (varian_read_parms(fname_parm, 1,
          VARIAN_PARMTYPE_FLOAT, parmstr, &parm_rfp) != 1)
      throw("failed to get %s from '%s'", parmstr, fname_parm);

    /* parse the nucleus name string. */
    snprintf(parmstr, N_BUF, "%cn%s", d > 0 ? 'd' : 't', d > 1 ? dstr : "");
    if (varian_read_parms(fname_parm, 1,
          VARIAN_PARMTYPE_STRING, parmstr, D->dims[d].nuc) != 1)
      throw("failed to get %s from '%s'", parmstr, fname_parm);

    /* store the read parameters in the datum structure. */
    D->dims[d].td = D->dims[d].sz = parm_np;
    D->dims[d].tdunif = parm_npmax;

    /* determine complexity based on dimension. */
    if (d == 0) {
      /* direct dimension: assume complex. */
      D->dims[d].cx = 1;
      D->dims[d].sz /= 2;
    }
    else {
      /* indirect dimension: parse the phase parameter. */
      snprintf(parmstr, N_BUF, "phase%s", d > 1 ? dstr : "");
      if (varian_read_parms(fname_parm, 1,
            VARIAN_PARMTYPE_INTS, parmstr, parm_phase) != 1)
        throw("failed to get %s from '%s'", parmstr, fname_parm);

      /* determine the phasing. */
      if (parm_phase[0] > 1) {
        /* the dimension is complex. */
        D->dims[d].cx = 1;
        D->dims[d].td *= 2;
        D->dims[d].tdunif *= 2;
      }
      else if (parm_phase[0] == 1 && parm_phase[1] == 3) {
        /* the dimension requires sign alternation. */
        D->dims[d].alt = 1;
      }
    }

    /* store the read spectral parameters in the datum structure. */
    D->dims[d].carrier = parm_sfrq;
    D->dims[d].width = parm_swh;
    D->dims[d].offset = parm_rfp;

    /* determine if the dimension is nonuniformly subsampled. */
    if (D->dims[d].tdunif && D->dims[d].td != D->dims[d].tdunif)
      D->dims[d].nus = 1;
  }

  /* check the ordering of the indirect dimensions. */
  if (D->nd >= 3) {
    /* parse the array parameter. */
    if (varian_read_parms(fname_parm, 1,
          VARIAN_PARMTYPE_STRING, "array", parmstr) != 1)
      throw("failed to parse array from '%s'", fname_parm);

    /* parse the array parameter. */
    arry = strsplit(parmstr, ",", &narry);

    /* check that the string was split. */
    if (!arry || !narry)
      throw("failed to split array string");

    /* allocate an array of dimension indices. */
    ord = (int*) calloc(D->nd, sizeof(int));
    if (!ord)
      throw("failed to allocate %u indices", D->nd);

    /* initialize the order array. */
    for (d = 0; d < D->nd; d++)
      ord[d] = -1;

    /* loop through the phase array. */
    for (i = 0; i < narry; i++) {
      /* determine the dimension index of the current phase. */
      for (d = 1; d < D->nd; d++) {
        /* build the dimension string and phase string. */
        snprintf(dstr, 8, "%u", d);
        snprintf(parmstr, N_BUF, "phase%s", d > 1 ? dstr : "");

        /* check if the current phase string matches the array value. */
        if (strcmp(parmstr, arry[i]) == 0) {
          /* yes. store it and break. */
          ord[i + 1] = d;
          break;
        }
      }
    }

    /* check if any values were not identified. */
    for (d = 0, i = 1; d < D->nd; d++)
      if (ord[d] == -1) i = 0;

    /* swap the ordering of the dimension parameters, but only if all
     * phase values were found in the array.
     */
    if (i && !datum_dims_reorder(D, ord))
      throw("failed to reorder dimensions");

    /* free the string and index arrays. */
    strvfree(arry, narry);
    free(ord);
  }

  /* build the schedule filename. */
  snprintf(fname_parm, n_fname, "%s/nls.hdr_3", dname);

  /* read the schedule file, if one exists. */
  if (bytes_fexist(fname_parm) && !datum_sched_read(D, fname_parm))
    throw("failed to read varian schedule");

  /* free the allocated filename strings. */
  D->fname = fname_data;
  free(fname_parm);

  /* store the datum type. */
  D->type = DATUM_TYPE_VARIAN;
  D->endian = BYTES_ENDIAN_AUTO;

  /* return success. */
  return 1;
}

