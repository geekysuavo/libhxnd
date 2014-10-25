
/* ndmath: A framework for n-dimensional hypercomplex calculations for NMR
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
#include <hxnd/varian.h>

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
      while (buf[strlen(buf) - 1] == '\n')
        buf[strlen(buf) - 1] = '\0';

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
int varian_read_hdr_file (const char *fname, unsigned int *endianness,
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
   * it's highly probable that the read bytes were big-endian and need
   * to be swapped.
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

    /* msb first. */
    *endianness = BYTES_ENDIAN_BIG;
  }
  else {
    /* lsb first. */
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
  unsigned int n, end, flt, nblk, szblk, offblk, offhead;
  struct varian_hdr_file fh;
  uint8_t *bytes;

  /* read the file header. */
  if (!varian_read_hdr_file(fname, &end, &fh))
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
  if (!bytes_toarray(bytes, n, end, fh.ebytes, flt, x))
    return 0;

  /* free the read byte data. */
  free(bytes);

  /* deinterlace the real and imaginary points into complex points. */
  if (!hx_array_deinterlace(x))
    return 0;

  /* return success. */
  return 1;
}

