
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

/* datum_dim_parms: array of dimension description structures that indicate
 * the names, types and offets of each datum dimension structure member .
 */
static const datum_dim_desc datum_dim_parms[] = {
  { "sz",      'u', offsetof(datum_dim, sz) },
  { "td",      'u', offsetof(datum_dim, td) },
  { "tdunif",  'u', offsetof(datum_dim, tdunif) },
  { "complex", 'u', offsetof(datum_dim, cx) },
  { "nus",     'u', offsetof(datum_dim, nus) },
  { "ft",      'u', offsetof(datum_dim, ft) },
  { "carrier", 'f', offsetof(datum_dim, carrier) },
  { "width",   'f', offsetof(datum_dim, width) },
  { "offset",  'f', offsetof(datum_dim, offset) },
  { "name",    's', offsetof(datum_dim, nuc) },
  { NULL, 'c', 0 }
};

/* datum_init(): initializes the elements of an NMR datum structure.
 * @D: pointer to the datum to initialize.
 */
void datum_init (datum *D) {
  /* initialize the input data filename. */
  D->fname = NULL;

  /* initialize the datum type. */
  D->type = DATUM_TYPE_UNDEFINED;

  /* initialize the dimension count and dimensions array. */
  D->dims = NULL;
  D->nd = 0;

  /* indicate that the array is not allocated. */
  D->array_alloc = 0;
}

/* datum_free(): frees all allocated innards of an NMR datum structure.
 * @D: pointer to the datum to free.
 */
void datum_free (datum *D) {
  /* free the input filename. */
  if (D->fname)
    free(D->fname);

  /* free the dimensions array. */
  if (D->dims)
    free(D->dims);

  /* free the array data. */
  datum_free_array(D);

  /* re-initialize the datum. */
  datum_init(D);
}

/* datum_guess_type(): attempts to reliably determine the type of file format
 * contained in a current file or directory.
 * @fname: the file or directory name to check.
 */
enum datum_type datum_guess_type (const char *fname) {
  /* check if the file is an hx-format file. */
  if (datum_check_magic(fname))
    return DATUM_TYPE_HXND;
  else
    traceback_clear();

  /* check if the file is a pipe-format file. */
  if (pipe_check_magic(fname))
    return DATUM_TYPE_PIPE;
  else
    traceback_clear();

  /* check if the file is a varian-format directory. */
  if (varian_check_dir(fname))
    return DATUM_TYPE_VARIAN;
  else
    traceback_clear();

  /* check if the file is a bruker-format directory. */
  if (bruker_check_dir(fname))
    return DATUM_TYPE_BRUKER;
  else
    traceback_clear();

  /* return no match. */
  return DATUM_TYPE_UNDEFINED;
}

/* datum_print(): prints the metadata associated with an acquired NMR datum.
 * @D: the datum to print data from.
 * @fname: the output filename.
 */
int datum_print (datum *D, const char *fname) {
  /* declare a few required variables:
   * @d: dimension loop counter.
   * @fh: the file handle used for writing.
   */
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

  /* print the first line. */
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

      /* print an ending if required. */
      if (d == D->array.k - 1)
        fprintf(fh, ")\n");
    }

    /* print a newline. */
    fprintf(fh, "\n");
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

  /* print the points count. */
  fprintf(fh, "Total:    ");
  for (d = 0; d < D->nd; d++)
    fprintf(fh, "%15u", D->dims[d].td);

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

/* datum_check_magic(): checks the first "magic" bytes of a file and
 * returns whether they match the NMR datum magic number.
 * @fname: the input filename.
 */
int datum_check_magic (const char *fname) {
  /* declare a few required variables:
   * @fh: input file handle.
   * @wd: first word of input file.
   */
  uint64_t wd;
  FILE *fh;

  /* open the input file. */
  fh = fopen(fname, "rb");

  /* check that the file was opened. */
  if (!fh)
    throw("failed to open '%s'", fname);

  /* read the first word. */
  if (!fread(&wd, sizeof(uint64_t), 1, fh))
    throw("failed to read magic number");

  /* close the input file. */
  fclose(fh);

  /* check the magic word, without swapping. */
  if (wd == NMR_DATUM_MAGIC)
    return 1;

  /* swap the word and check again. */
  bytes_swap_u64(&wd);
  if (wd == NMR_DATUM_MAGIC)
    return 1;

  /* no match. */
  return 0;
}

/* define the number of (u64) members in the header and dimension sections of
 * binary datum files.
 */
#define NMR_DATUM_FWRITE_SZ_HDR  4
#define NMR_DATUM_FWRITE_SZ_DIM  8

/* define bit field positions to store status flags in binary datum files.
 */
#define NMR_DATUM_S_COMPLEX  0x0000000000000001
#define NMR_DATUM_S_NUS      0x0000000000000002
#define NMR_DATUM_S_FFT      0x0000000000000004

/* datum_fwrite(): writes an NMR datum structure to an opened file stream.
 * @D: pointer to the source structure.
 * @fh: the output file stream.
 */
int datum_fwrite (datum *D, FILE *fh) {
  /* declare a few required variables:
   * @i: buffer index.
   * @d: dimension loop counter.
   * @n_buf: current buffer word size.
   * @buf: output header/dimension buffer.
   * @status: status word value.
   */
  unsigned int i, d, n_buf;
  uint64_t *buf, status;

  /* allocate the header buffer. */
  n_buf = NMR_DATUM_FWRITE_SZ_HDR;
  buf = (uint64_t*) calloc(n_buf, sizeof(uint64_t));

  /* check that allocation succeeded. */
  if (!buf)
    throw("failed to allocate header buffer");

  /* build the file header. */
  i = 0;
  buf[i++] = (uint64_t) NMR_DATUM_MAGIC;
  buf[i++] = (uint64_t) D->endian;
  buf[i++] = (uint64_t) DATUM_TYPE_HXND;
  buf[i++] = (uint64_t) D->nd;

  /* write the file header. */
  if (fwrite(buf, sizeof(uint64_t), n_buf, fh) != n_buf)
    throw("failed to write header");

  /* reallocate the buffer to hold dimension information. */
  n_buf = NMR_DATUM_FWRITE_SZ_DIM;
  buf = (uint64_t*) realloc(buf, n_buf * sizeof(uint64_t));

  /* check that reallocation succeeded. */
  if (!buf)
    throw("failed to allocate dimension buffer");

  /* loop over the dimensions in the datum structure. */
  for (d = 0; d < D->nd; d++) {
    /* zero the buffer memory. */
    memset(buf, 0, n_buf * sizeof(uint64_t));

    /* store the size parameters. */
    i = 0;
    buf[i++] = (uint64_t) D->dims[d].sz;
    buf[i++] = (uint64_t) D->dims[d].td;
    buf[i++] = (uint64_t) D->dims[d].tdunif;

    /* build the status word. */
    status = 0;
    status |= (D->dims[d].cx ?  NMR_DATUM_S_COMPLEX : 0);
    status |= (D->dims[d].nus ? NMR_DATUM_S_NUS : 0);
    status |= (D->dims[d].ft ?  NMR_DATUM_S_FFT : 0);

    /* store the status word. */
    buf[i++] = status;

    /* store the spectral parameters. */
    buf[i++] = bytes_real_to_u64(D->dims[d].carrier);
    buf[i++] = bytes_real_to_u64(D->dims[d].width);
    buf[i++] = bytes_real_to_u64(D->dims[d].offset);

    /* store the nucleus string. */
    memcpy(buf + (i++), D->dims[d].nuc, sizeof(uint64_t));

    /* write the dimension buffer. */
    if (fwrite(buf, sizeof(uint64_t), n_buf, fh) != n_buf)
      throw("failed to write dimension %u", d);
  }

  /* free the allocated buffer. */
  free(buf);

  /* write the core array content to the end of the file stream. */
  if (D->array_alloc && !hx_array_fwrite(&D->array, fh))
    throw("failed to write core array");

  /* return success. */
  return 1;
}

/* datum_fread(): reads an NMR datum structure from an opened file stream.
 * @D: pointer to the destination structure.
 * @fh: the input file stream.
 * @read_array: whether to read the array data.
 */
int datum_fread (datum *D, FILE *fh, int read_array) {
  /* declare a few required variables:
   * @i: buffer index.
   * @d: dimension loop counter.
   * @n_buf: current buffer word size.
   * @buf: output header/dimension buffer.
   * @status: status word value.
   */
  unsigned int i, d, n_buf, swapping;
  uint64_t *buf, status;

  /* allocate the header buffer. */
  n_buf = NMR_DATUM_FWRITE_SZ_HDR;
  buf = (uint64_t*) calloc(n_buf, sizeof(uint64_t));

  /* check that allocation succeeded. */
  if (!buf)
    throw("failed to allocate header buffer");

  /* read the file header. */
  if (fread(buf, sizeof(uint64_t), n_buf, fh) != n_buf)
    throw("failed to read header");

  /* check the first word in the header. if it does not match, then
   * byte-swap the data.
   */
  if (buf[0] != NMR_DATUM_MAGIC) {
    /* no match. swap the bytes of each word. */
    bytes_swap_general((uint8_t*) buf,
                       n_buf * sizeof(uint64_t),
                       sizeof(uint64_t));

    /* now check the magic word. */
    if (buf[0] != NMR_DATUM_MAGIC)
      throw("invalid magic word 0x%08x", buf[0]);

    /* set the swapping flag. */
    swapping = 1;
  }
  else {
    /* match. no swaps needed. */
    swapping = 0;
  }

  /* unpack the file header. */
  i = 1;
  D->endian = (enum byteorder) buf[i++];
  D->type = (enum datum_type) buf[i++];
  D->nd = (unsigned int) buf[i++];

  /* allocate the dimension array. */
  D->dims = (datum_dim*) calloc(D->nd, sizeof(datum_dim));

  /* check that allocation succeeded. */
  if (D->dims == NULL)
    throw("failed to allocate %u dimensions", D->nd);

  /* reallocate the buffer to hold dimension information. */
  n_buf = NMR_DATUM_FWRITE_SZ_DIM;
  buf = (uint64_t*) realloc(buf, n_buf * sizeof(uint64_t));

  /* check that reallocation succeeded. */
  if (!buf)
    throw("failed to allocate dimension buffer");

  /* loop over the dimensions to be read in. */
  for (d = 0; d < D->nd; d++) {
    /* read the dimension data. */
    if (fread(buf, sizeof(uint64_t), n_buf, fh) != n_buf)
      throw("failed to read dimension %u", d);

    /* swap the bytes, if required. */
    if (swapping)
      bytes_swap_general((uint8_t*) buf,
                         n_buf * sizeof(uint64_t),
                         sizeof(uint64_t));

    /* unpack the dimension size parameters. */
    i = 0;
    D->dims[d].sz = (unsigned int) buf[i++];
    D->dims[d].td = (unsigned int) buf[i++];
    D->dims[d].tdunif = (unsigned int) buf[i++];

    /* unpack the status word. */
    status = buf[i++];
    D->dims[d].cx = (status & NMR_DATUM_S_COMPLEX ? 1 : 0);
    D->dims[d].nus = (status & NMR_DATUM_S_NUS ? 1 : 0);
    D->dims[d].ft = (status & NMR_DATUM_S_FFT ? 1 : 0);

    /* unpack the spectral parameters. */
    D->dims[d].carrier = bytes_u64_to_real(buf[i++]);
    D->dims[d].width = bytes_u64_to_real(buf[i++]);
    D->dims[d].offset = bytes_u64_to_real(buf[i++]);

    /* unpack the nucleus string. */
    memcpy(D->dims[d].nuc, buf + (i++), sizeof(uint64_t));
  }

  /* free the allocated buffer. */
  free(buf);

  /* check if the core array content should be loaded. */
  if (read_array) {
    /* read the core array content from the end of the file stream. */
    if (!hx_array_fread(&D->array, fh))
      throw("failed to read core array");

    /* read succeeded. set the array allocation status flag. */
    D->array_alloc = 1;
  }

  /* return success. */
  return 1;
}

/* datum_save(): saves the contents of an acquired NMR datum to a file,
 * or standard output if a NULL filename was passed.
 * @D: the datum to save data from.
 * @fname: the output filename.
 */
int datum_save (datum *D, const char *fname) {
  /* declare a required variable. */
  FILE *fh;

  /* open the output file. */
  if (fname)
    fh = fopen(fname, "wb");
  else
    fh = stdout;

  /* check that the file was opened. */
  if (!fh)
    throw("failed to open '%s'", fname);

  /* write the datum to the file. */
  if (!datum_fwrite(D, fh))
    throw("failed to write '%s'", fname);

  /* close the output file. */
  if (fname)
    fclose(fh);

  /* return success. */
  return 1;
}

/* datum_load(): loads the contents of an acquired NMR datum from a file.
 * @D: the datum to load data from.
 * @fname: the input filename.
 * @load_array: whether to load the array data.
 */
int datum_load (datum *D, const char *fname, int load_array) {
  /* declare a required variable. */
  FILE *fh;

  /* open the input file. */
  if (fname)
    fh = fopen(fname, "rb");
  else
    fh = stdin;

  /* check that the file was opened. */
  if (!fh)
    throw("failed to open '%s'", fname);

  /* read the datum from the file. */
  if (!datum_fread(D, fh, load_array))
    throw("failed to read '%s'", fname);

  /* because we're here, we know the read succeeded, so check if the data
   * was loaded from a named file.
   */
  if (fname) {
    /* store the input filename. */
    D->fname = (char*) malloc((strlen(fname) + 1) * sizeof(char));
    if (D->fname)
      strcpy(D->fname, fname);

    /* close the input file. */
    fclose(fh);
  }

  /* return success. */
  return 1;
}

/* datum_get_dim_parameter(): retrieves a dimension parameter from a datum
 * structure, by name and dimension index.
 * @D: pointer to the datum struct.
 * @name: parameter name to get.
 * @d: dimension of parameter.
 * @parm: pointer to the result.
 */
int datum_get_dim_parameter (datum *D, const char *name, unsigned int d,
                             void *parm) {
  /* declare pointer locations for each data type. */
  unsigned int *uptr;
  char *sptr;
  real *fptr;
  int *dptr;

  /* declare required variables. */
  unsigned int i;
  uint8_t *base;

  /* check that the dimension index is within bounds. */
  if (d >= D->nd)
    throw("dimension index %u out of bound %u", d, D->nd);

  /* compute the base address of the dimension structure of interest. */
  base = (uint8_t*) &D->dims[d];

  /* loop over the dimension descriptors. */
  for (i = 0; datum_dim_parms[i].name; i++) {
    /* check if the current descriptor matches the requested name. */
    if (strcmp(datum_dim_parms[i].name, name) == 0) {
      /* act based on the parameter type. */
      switch (datum_dim_parms[i].type) {
        /* signed int */
        case 'd':
          dptr = (int*) (base + datum_dim_parms[i].off);
          *((int*) parm) = *dptr;
          return 1;

        /* unsigned int */
        case 'u':
          uptr = (unsigned int*) (base + datum_dim_parms[i].off);
          *((unsigned int*) parm) = *uptr;
          return 1;

        /* real */
        case 'f':
          fptr = (real*) (base + datum_dim_parms[i].off);
          *((real*) parm) = *fptr;
          return 1;

        /* string */
        case 's':
          sptr = (char*) (base + datum_dim_parms[i].off);
          memcpy(parm, sptr, 8);
          return 1;

        /* unknown */
        default:
          throw("invalid parameter type '%c'", datum_dim_parms[i].type);
      }
    }
  }

  /* no match was found. throw an exception. */
  throw("invalid dimension parameter name '%s'", name);
}

/* datum_set_dim_parameter(): stores a dimension parameter into a datum
 * structure, by name and dimension index.
 * @D: pointer to the datum struct.
 * @name: parameter name to set.
 * @d: dimension of parameter.
 * @parm: the new value, as a string.
 */
int datum_set_dim_parameter (datum *D, const char *name, unsigned int d,
                             const char *parm) {
  /* declare pointer locations for each data type. */
  unsigned int *uptr;
  char *sptr;
  real *fptr;
  int *dptr;

  /* declare required variables. */
  unsigned int i;
  uint8_t *base;

  /* check that the dimension index is within bounds. */
  if (d >= D->nd)
    throw("dimension index %u out of bound %u", d, D->nd);

  /* compute the base address of the dimension structure of interest. */
  base = (uint8_t*) &D->dims[d];

  /* loop over the dimension descriptors. */
  for (i = 0; datum_dim_parms[i].name; i++) {
    /* check if the current descriptor matches the requested name. */
    if (strcmp(datum_dim_parms[i].name, name) == 0) {
      /* act based on the parameter type. */
      switch (datum_dim_parms[i].type) {
        /* signed int */
        case 'd':
          dptr = (int*) (base + datum_dim_parms[i].off);
          *dptr = atol(parm);
          return 1;

        /* unsigned int */
        case 'u':
          uptr = (unsigned int*) (base + datum_dim_parms[i].off);
          *uptr = (unsigned int) atol(parm);
          return 1;

        /* real */
        case 'f':
          fptr = (real*) (base + datum_dim_parms[i].off);
          *fptr = (real) atof(parm);
          return 1;

        /* string */
        case 's':
          sptr = (char*) (base + datum_dim_parms[i].off);
          memcpy(sptr, parm, 8);
          return 1;

        /* unknown */
        default:
          throw("invalid parameter type '%c'", datum_dim_parms[i].type);
      }
    }
  }

  /* no match was found. throw an exception. */
  throw("invalid dimension parameter name '%s'", name);
}

/* datum_reorder_dims(): reorders the dimensions inside a datum structure
 * according to the initial ordering @order. the target ordering will be
 * achieved by sorting the values in @order until they are in increasing
 * order.
 *
 * @D: pointer to the target datum structure.
 * @order: initial dimension ordering.
 */
int datum_reorder_dims (datum *D, int *order) {
  /* declare required variables:
   * @i: insertion sort outer loop counter.
   * @j: insertion sort inner loop counter.
   * @dim: a swap location for dimension information.
   */
  unsigned int i, j;
  datum_dim dim;
  int *ord, swp;

  /* allocate a temporary array of indices. */
  ord = (int*) malloc(D->nd * sizeof(int));
  if (!ord)
    throw("failed to allocate %u indices", D->nd);

  /* copy the ordering into the temporary array. */
  memcpy(ord, order, D->nd * sizeof(int));

  /* loop over the array of dimensions. */
  for (i = 1; i < D->nd; i++) {
    /* set the initial inner loop index. */
    j = i;

    /* loop over the unsorted dimensions. */
    while (j > 0 && ord[j - 1] > ord[j]) {
      /* swap the (j-1) and (j) elements. */
      dim = D->dims[j];
      D->dims[j] = D->dims[j - 1];
      D->dims[j - 1] = dim;

      /* swap the ordering of the indices. */
      swp = ord[j];
      ord[j] = ord[j - 1];
      ord[j - 1] = swp;

      /* decrement the inner loop counter. */
      j--;
    }
  }

  /* free the ordering array. */
  free(ord);

  /* return success. */
  return 1;
}

/* datum_refactor_array(): repacks, infills and deinterlaces the core array
 * structure of an NMR datum until it's dimensionality and complexity agree
 * with the dimension parameter values.
 * @D: pointer to the datum to manipulate.
 */
int datum_refactor_array (datum *D) {
  /* declare a few required variables:
   * @d: dimension loop counter.
   */
  unsigned int d;

  /* check that the array has been allocated. */
  if (!D->array_alloc)
    throw("array is unallocated");

  /* loop over the acquisition dimensions to refactor the nD array. */
  for (d = 0; d < D->nd; d++) {
    /* repack indirect dimensions in the array. */
    if (d > 0 && !hx_array_repack(&D->array, D->dims[d - 1].sz))
      throw("failed to repack array dimension %d", d);

    /* check if the current dimension is complex. */
    if (D->dims[d].cx) {
      /* de-interlace this dimension. */
      if (!hx_array_deinterlace(&D->array))
        throw("failed to deinterlace complex dimension %d", d);
    }
    else {
      /* increment the dimensionality without de-interlacing. */
      if (!hx_array_resize(&D->array,
            D->array.d + 1, D->array.k, D->array.sz))
        throw("failed to complex-promote dimension %d", d);
    }
  }

  /* return success. */
  return 1;
}

/* datum_read_array(): reads and refactors the array data of a datum structure
 * that has been initialized with a correct set of parameters.
 * @D: pointer to the datum to manipulate.
 */
int datum_read_array (datum *D) {
  /* declare a few required variables. */
  unsigned int d, nblk, szblk, offset;
  FILE *fh;

  /* check if the array has been allocated. */
  if (D->array_alloc)
    throw("array is already allocated");

  /* check that the filename is non-null and non-empty. */
  if (D->fname == NULL || strlen(D->fname) == 0)
    throw("filename is invalid");

  /* load based on the type of data. */
  if (D->type == DATUM_TYPE_BRUKER) {
    /* determine the data block size. */
    szblk = 4 * D->dims[0].td;

    /* determine the data block count. */
    for (d = 1, nblk = 1; d < D->nd; d++)
      nblk *= D->dims[d].td;

    /* check if the blocks are 1.0 KiB-aligned. */
    if (szblk % 1024 == 0) {
      /* yes. use a (faster) single read, because no gaps exist. */
      szblk *= nblk;
      nblk = 1;
    }

    /* load the raw data from the fid/ser file. */
    if (!bruker_read(D->fname, D->endian, nblk, szblk, &D->array))
      throw("failed to read bruker data from '%s'", D->fname);
  }
  else if (D->type == DATUM_TYPE_VARIAN) {
    /* load the raw data from the fid file. */
    if (!varian_read(D->fname, &D->array))
      throw("failed to read varian data from '%s'", D->fname);
  }
  else if (D->type == DATUM_TYPE_PIPE) {
    /* determine the data byte count. */
    for (d = 0, szblk = sizeof(float); d < D->nd; d++)
      szblk *= (D->dims[d].cx ? 2 : 1) * D->dims[d].sz;

    /* load the raw data from the pipe file. */
    if (!pipe_read(D->fname, szblk, &D->array))
      throw("failed to read pipe data from '%s'", D->fname);

    /* interlace the complex points, if required. */
    if (D->dims[0].cx && !pipe_interlace(&D->array, D->dims[0].sz))
      throw("failed to interlace complex traces");
  }
  else if (D->type == DATUM_TYPE_HXND) {
    /* compute the offset where the array begins. */
    offset = NMR_DATUM_FWRITE_SZ_HDR;
    offset += D->nd * NMR_DATUM_FWRITE_SZ_DIM;
    offset *= sizeof(uint64_t);

    /* open the input file. */
    fh = fopen(D->fname, "rb");

    /* check that the file was opened. */
    if (!fh)
      throw("failed to open '%s'", D->fname);

    /* seek past the file header information. */
    if (fseek(fh, offset, SEEK_SET))
      throw("failed to seek to array in '%s'", D->fname);

    /* read the array data from the hx-format file. */
    if (!hx_array_fread(&D->array, fh))
      throw("failed to read hx-format data from '%s'", D->fname);

    /* close the input file. */
    fclose(fh);
  }
  else
    throw("unsupported data type %d", D->type);

  /* indicate that the array has been allocated. */
  D->array_alloc = 1;

  /* refactor the core array, but only if it's *not* already in
   * the native hxnd format.
   */
  if (D->type != DATUM_TYPE_HXND && !datum_refactor_array(D))
    throw("failed to refactor array");

  /* return success. */
  return 1;
}

/* datum_free_array(): frees an allocated array structure from an NMR datum.
 * @D: pointer to the datum to manipulate.
 */
int datum_free_array (datum *D) {
  /* check if the array is allocated. */
  if (!D->array_alloc)
    return 1;

  /* free the data array. */
  hx_array_free(&D->array);

  /* indicate that the array has been de-allocated. */
  D->array_alloc = 0;

  /* return success. */
  return 1;
}

/* datum_fill(): parses acquisition parameters into an NMR datum structure.
 * the 'type' member of the datum structure must be initialized to the type
 * of parms/data to parse in.
 * @D: pointer to the datum struct to fill.
 * @fname: input file or directory name.
 */
int datum_fill (datum *D, const char *fname) {
  /* act based on the datum type. */
  switch (D->type) {
    /* undefined. */
    case DATUM_TYPE_UNDEFINED:
      /* throw an error for unknown data types. */
      throw("undefined data type");

    /* bruker. */
    case DATUM_TYPE_BRUKER:
      /* attempt to fill from bruker-format files. */
      if (!bruker_fill_datum(fname, D))
        throw("failed to parse bruker parameters");

      /* break out. */
      break;

    /* varian. */
    case DATUM_TYPE_VARIAN:
      /* attempt to fill from varian-format files. */
      if (!varian_fill_datum(fname, D))
        throw("failed to parse varian parameters");

      /* break out. */
      break;

    /* pipe. */
    case DATUM_TYPE_PIPE:
      /* attempt to fill from a pipe-format file. */
      if (!pipe_fill_datum(fname, D))
        throw("failed to parse pipe parameters");

      /* break out. */
      break;

    /* hxnd (native). */
    case DATUM_TYPE_HXND:
      /* just load the datum file. */
      if (!datum_load(D, fname, 0))
        throw("failed to parse hx parameters");

      /* break out. */
      break;

    /* other. */
    default:
      throw("unknown data type %d", D->type);
  }

  /* return success. */
  return 1;
}

