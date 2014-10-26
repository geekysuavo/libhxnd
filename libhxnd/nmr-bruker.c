
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
   * @bufkey: parsed key string from each buffer line.
   * @bufval: parsed value string from each buffer line.
   * @peq: pointer inside @buf to the equals sign.
   * @ncpy: number of chars to copy for key/value pairs.
   * @fh: the input file handle.
   */
  char buf[N_BUF], bufkey[N_BUF], bufval[N_BUF], *peq;
  unsigned int ncpy;
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

      /* check if the current line can contain a key/value pair. */
      if (strlen(buf) <= 3 || strncmp(buf, "##$", 3))
        continue;

      /* locate the equals sign in the current line. */
      peq = strstr(buf, "= ");

      /* check if the line contains an equals sign. */
      if (!peq)
        continue;

      /* copy the key substring from the current line. */
      ncpy = peq - buf - 3;
      strncpy(bufkey, buf + 3, ncpy);
      bufkey[ncpy] = '\0';

      /* copy the value substring from the current line. */
      ncpy = strlen(buf) - ncpy - 6;
      strncpy(bufval, peq + 2, ncpy);
      bufval[ncpy] = '\0';

      /* loop through the parameters to search for. */
      for (i = 0; i < n; i++) {
        /* check if the current parameter key is a match. */
        if (strcmp(keys[i], bufkey) == 0) {
          /* yep! increment the number of identified parameters. */
          nid++;

          /* act based on the parameter type. */
          switch (typs[i]) {
            /* integer. */
            case BRUKER_PARMTYPE_INT:
              *((int*) vals[i]) = atol(bufval);
              break;

            /* float */
            case BRUKER_PARMTYPE_FLOAT:
              *((float*) vals[i]) = atof(bufval);
              break;

            /* string */
            case BRUKER_PARMTYPE_STRING:
              ncpy = strlen(bufval) - 2;
              strncpy(vals[i], bufval + 1, ncpy);
              ((char*) vals[i])[ncpy] = '\0';
              break;

            /* other. */
            default:
              break;
          }
        }
      }
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

/* bruker_read(): reads a bruker data file into a complex linear array.
 * @fname: the input data filename.
 * @endianness: the data byte ordering.
 * @nblk: the number of blocks/fids.
 * @szblk: the size of block/fid.
 * @x: the output array.
 */
int bruker_read (const char *fname, unsigned int endianness,
                 unsigned int nblk, unsigned int szblk,
                 hx_array *x) {
  /* declare a few required variables.
   */
  unsigned int n;
  uint8_t *bytes;

  /* read the data bytes in from the file. */
  bytes = bytes_read_bruker(fname, nblk, szblk, &n);

  /* check that the bytes were read successfully. */
  if (!bytes)
    return 0;

  /* build a real linear array from the byte data. */
  if (!bytes_toarray(bytes, n, endianness, 4, 0, x))
    return 0;

  /* free the read byte data. */
  free(bytes);

  /* deinterlace the real and imaginary points into complex points. */
  if (!hx_array_deinterlace(x))
    return 0;

  /* return success. */
  return 1;
}

