
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
    return 0;

  /* initialize the variable arguments list. */
  va_start(vl, n);

  /* allocate memory to store the search options for each parameter. */
  keys = (char**) malloc(n * sizeof(char*));
  vals = (void**) malloc(n * sizeof(void*));
  typs = (char*) malloc(n * sizeof(char));

  /* check that the memory was allocated successfully. */
  if (!keys || !vals || !typs)
    return 0;

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
      return 0;

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
            while (buf[strlen(buf) - 1] == '\n')
              buf[strlen(buf) - 1] = '\0';

            /* split the next line by whitespace. */
            fields = strsplit(buf, " ", &nfields);

            /* act based on the parameter type. */
            switch (typs[i]) {
              /* integer. */
              case VARIAN_PARMTYPE_INT:
                *((int*) vals[i]) = atol(fields[1]);
                break;

              /* float */
              case VARIAN_PARMTYPE_FLOAT:
                *((float*) vals[i]) = atof(fields[1]);
                break;

              /* string */
              case VARIAN_PARMTYPE_STRING:
                ncpy = strlen(fields[1]) - 2;
                strncpy(vals[i], fields[1] + 1, ncpy);
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

/* varian_read_hdr_file(): reads a varian data file header into a struct.
 * @fname: the data filename.
 * @endianness: byte order result pointer.
 * @hdr: header result pointer.
 */
int varian_read_hdr_file (const char *fname, enum byteorder *endianness,
                          struct varian_hdr_file *hdr) {
  /* declare a few required variables. */
  uint8_t *bytes;

  /* read in the file header bytes. */
  bytes = bytes_read_block(fname, 0, sizeof(struct varian_hdr_file));

  /* check that the bytes were read successfully. */
  if (!bytes)
    return 0;

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
    if (bytes_native(BYTES_ENDIAN_BIG))
      *endianness = BYTES_ENDIAN_LITTLE;
    else
      *endianness = BYTES_ENDIAN_BIG;
  }
  else {
    /* same endianness. */
    if (bytes_native(BYTES_ENDIAN_BIG))
      *endianness = BYTES_ENDIAN_BIG;
    else
      *endianness = BYTES_ENDIAN_LITTLE;
  }

  /* return success. */
  return 1;
}

/* varian_read(): reads a varian data file into a complex linear array.
 * @fname: the input data filename.
 * @x: the output array.
 */
int varian_read (const char *fname, hx_array *x) {
  /* declare a few required variables.
   */
  unsigned int n, flt, nblk, szblk, offblk, offhead;
  struct varian_hdr_file fh;
  enum byteorder endianness;
  uint8_t *bytes;

  /* read the file header. */
  if (!varian_read_hdr_file(fname, &endianness, &fh))
    return 0;

  /* compute block and header sizes for data loading. */
  nblk = fh.nblocks;
  szblk = fh.ntraces * fh.tbytes;
  offblk = fh.bbytes - szblk;
  offhead = sizeof(struct varian_hdr_file);

  /* determine whether the data is floating point or not. */
  flt = (fh.status & VARIAN_HDR_S_FLOAT ? 1 : 0);

  /* read the data bytes in from the file. */
  bytes = bytes_read_varian(fname, nblk, szblk, offblk, offhead, &n);

  /* check that the bytes were read successfully. */
  if (!bytes)
    return 0;

  /* build a real linear array from the byte data. */
  if (!bytes_toarray(bytes, n, endianness, fh.ebytes, flt, x))
    return 0;

  /* free the read byte data. */
  free(bytes);

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

/* varian_datum(): completely loads varian raw data into an NMR datum
 * structure.
 * @dname: the input directory name.
 * @D: pointer to the datum struct to fill.
 */
int varian_datum (const char *dname, datum *D) {
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
   * @dstr: dimension loop counter (string).
   */
  unsigned int d;
  char dstr[8];

  /* declare variables for parsing dimension parameters:
   * @parmstr: parameter string.
   */
  real parm_swh, parm_sfrq, parm_rfp;
  char parmstr[N_BUF];
  int parm_np;

  /* allocate memory for the fid and procpar filenames. */
  n_fname = strlen(dname) + 12;
  fname_data = (char*) malloc(n_fname * sizeof(char));
  fname_parm = (char*) malloc(n_fname * sizeof(char));

  /* check that the filename strings were allocated. */
  if (!fname_data || !fname_parm)
    return 0;

  /* build the procpar filename. */
  snprintf(fname_parm, n_fname, "%s/procpar", dname);

  /* build the fid filename. */
  snprintf(fname_data, n_fname, "%s/fid", dname);

  /* algorithmically guess the dimensionality of the data. */
  D->nd = varian_count_dims(fname_parm);

  /* check the dimensionality. */
  if (D->nd < 1)
    return 0;

  /* allocate the dimension parameter array. */
  D->dims = (datum_dim*) calloc(D->nd, sizeof(datum_dim));

  /* check that the dimension parameter array was allocated. */
  if (D->dims == NULL)
    return 0;

  /* loop over the acquisition dimensions. */
  for (d = 0; d < D->nd; d++) {
    /* build a string version of the dimension index. */
    snprintf(dstr, 8, "%u", d);

    /* parse the number of points. */
    snprintf(parmstr, N_BUF, "n%c%s", d > 0 ? 'i' : 'p', d > 1 ? dstr : "");
    if (varian_read_parms(fname_parm, 1,
          VARIAN_PARMTYPE_INT, parmstr, &parm_np) != 1)
      return 0;

    /* parse the carrier base frequency. */
    snprintf(parmstr, N_BUF, "%cfrq%s", d > 0 ? 'd' : 's', d > 1 ? dstr : "");
    if (varian_read_parms(fname_parm, 1,
          VARIAN_PARMTYPE_FLOAT, parmstr, &parm_sfrq) != 1)
      return 0;

    /* parse the spectral width. */
    snprintf(parmstr, N_BUF, "sw%s", d > 0 ? dstr : "");
    if (varian_read_parms(fname_parm, 1,
          VARIAN_PARMTYPE_FLOAT, parmstr, &parm_swh) != 1)
      return 0;

    /* parse the spectral offset. */
    snprintf(parmstr, N_BUF, "rfp%s", d > 0 ? dstr : "");
    if (varian_read_parms(fname_parm, 1,
          VARIAN_PARMTYPE_FLOAT, parmstr, &parm_rfp) != 1)
      return 0;

    /* parse the nucleus name string. */
    snprintf(parmstr, N_BUF, "%cn%s", d > 0 ? 'd' : 't', d > 1 ? dstr : "");
    if (varian_read_parms(fname_parm, 1,
          VARIAN_PARMTYPE_STRING, parmstr, D->dims[d].nuc) != 1)
      return 0;

    /* FIXME: correct procpar loading in varian_datum() */

    /* store the read parameters in the datum structure. */
    D->dims[d].td = D->dims[d].sz = parm_np;

    /* store the read spectral parameters in the datum structure. */
    D->dims[d].carrier = parm_sfrq;
    D->dims[d].width = parm_swh;
    D->dims[d].offset = parm_rfp;
  }

  /* FIXME: implement dimension re-ordering in varian_datum() */

  /* load the raw data from the fid file. */
  if (!varian_read(fname_data, &D->array))
    return 0;

  /* loop over the acquisition dimensions to refactor the nD array. */
  for (d = 0; d < D->nd; d++) {
    /* FIXME: implement array refactoring in varian_datum() */
  }

  /* free the allocated filename strings. */
  free(fname_data);
  free(fname_parm);

  /* return success. */
  return 1;
}

