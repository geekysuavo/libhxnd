
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

/* datum_type_def: datum type definition structure for holding all available
 * datum type names and values.
 */
struct datum_type_def {
  /* @name: the string name of the datum type.
   * @type: the value of the datum type.
   */
  const char *name;
  enum datum_type type;
};

/* datum_dim_desc: structure that defines a mapping between dimension
 * parameter names, their locations within the datum_dim structure,
 * and their types in the structure.
 */
typedef struct {
  /* @name: parameter name string.
   * @type: parameter type character.
   * @off: struct member offset.
   */
  const char *name;
  char type;
  size_t off;
}
datum_dim_desc;

/* datum_types: local table of all available datum type names and values.
 */
static const struct datum_type_def datum_types[] = {
  { "bruker", DATUM_TYPE_BRUKER },
  { "varian", DATUM_TYPE_VARIAN },
  { "pipe",   DATUM_TYPE_PIPE },
  { "hx",     DATUM_TYPE_HXND },
  { "text",   DATUM_TYPE_TEXT },
  { NULL, DATUM_TYPE_UNDEFINED }
};

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

  /* initialize the sampling schedule. */
  D->sched = NULL;
  D->d_sched = 0;
  D->n_sched = 0;

  /* initialize the group delay value. */
  D->grpdelay = 0.0;

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

  /* free the sampling schedule. */
  if (D->sched)
    free(D->sched);

  /* free the array data. */
  datum_free_array(D);

  /* re-initialize the datum. */
  datum_init(D);
}

/* datum_lookup_type(): returns the enumerated type based on its string
 * representation.
 * @name: the datum type string.
 */
enum datum_type datum_lookup_type (const char *name) {
  /* declare a required variable. */
  unsigned int i;

  /* loop over all supported datum types. */
  for (i = 0; datum_types[i].name; i++) {
    /* break if the datum name matches. */
    if (strcmp(name, datum_types[i].name) == 0)
      break;
  }

  /* return the identified type value. */
  return datum_types[i].type;
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

/* datum_fwrite_formatted(): writes an NMR datum structure an an opened file
 * stream in the specified output format.
 * @D: pointer to the source structure.
 * @fh: the output file stream.
 * @fmt: the output format.
 */
int datum_fwrite_formatted (datum *D, FILE *fh, enum datum_type fmt) {
  /* determine which write routine to utilize. */
  switch (fmt) {
    /* hx native format. */
    case DATUM_TYPE_HXND:
      return datum_fwrite(D, fh);

    /* text format. */
    case DATUM_TYPE_TEXT:
      return datum_fwrite_text(D, fh);

    /* default: unsupported. */
    default:
      throw("unsupported output format %d", fmt);
  }

  /* return success. */
  return 1;
}

/* datum_fwrite_text(): writes an NMR datum structure in text-format to
 * an opened file stream.
 * @D: pointer to the source structure.
 * @fh: the output file stream.
 */
int datum_fwrite_text (datum *D, FILE *fh) {
  /* declare required variables:
   * @d: dimension loop counter.
   */
  int i, *arr, idx;
  unsigned int d;

  /* print heading information for each dimension. */
  for (d = 0; d < D->nd; d++) {
    /* print the dimension index and name. */
    fprintf(fh, "# Axis %2u ('%s'):\n", d + 1, D->dims[d].nuc);

    /* print the dimension sizes. */
    fprintf(fh, "# Points:   %15u\n", D->dims[d].sz);
    fprintf(fh, "# Total:    %15u\n", D->dims[d].td);

    /* print the spectral parameters. */
    fprintf(fh, "# Obs (MHz):%15.3f\n", D->dims[d].carrier);
    fprintf(fh, "# SW (Hz):  %15.3f\n", D->dims[d].width);
    fprintf(fh, "# Off (Hz): %15.3f\n", D->dims[d].offset);

    /* print a spacer. */
    fprintf(fh, "#\n");
  }

  /* allocate an index array for iteration. */
  arr = hx_array_index_alloc(D->array.k);

  /* check that allocation was successful. */
  if (!arr)
    throw("failed to allocate %d indices", D->array.k);

  /* iterate over the points in the core array. */
  idx = 0;
  do {
    /* print the indices. */
    for (i = 0; i < D->array.k; i++)
      fprintf(fh, "%6d ", arr[i]);

    /* print the coefficients. */
    for (i = 0; i < D->array.n; i++)
      fprintf(fh, "%18.8e ", D->array.x[i + D->array.n * idx]);

    /* print a newline. */
    fprintf(fh, "\n");

    /* increment the linear index. */
    idx++;
  } while (hx_array_index_incr(D->array.k, D->array.sz, arr));

  /* free the index array. */
  free(arr);

  /* return success. */
  return 1;
}

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

  /* FIXME: implement schedule writing in datum_fwrite() */

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
    bytes_swap((uint8_t*) buf, n_buf, sizeof(uint64_t));

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
      bytes_swap((uint8_t*) buf, n_buf, sizeof(uint64_t));

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

  /* FIXME: implement schedule reading in datum_fread() */

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

/* datum_infill_array(): internally zero-fills the core array structure of an
 * NMR datum to reflect the final Nyquist grid. if no datum dimensions are
 * nonuniformly sampled, this routine has no effect.
 * @D: pointer to the datum to manipulate.
 *
 * NOTE: this routine requires that D->dims[0].sz matches D->array.sz[0].
 */
int datum_infill_array (datum *D) {
  /* declare a few required variables:
   * @i: loop counter for sampling schedule rows.
   * @j: loop counter for sampling schedule columns.
   * @idxi: input array linear index.
   * @idxo: output array linear index.
   * @arr: output array unpacked indices.
   * @nnus: number of nonuniform dimensions.
   * @ncpy: number of array points per sampled trace.
   * @nbytes: number of bytes per sampled trace.
   * @sznew: output array point count.
   * @tdnew: output array unpacked sizes.
   * @d: datum dimension loop counter.
   * @anew: output array structure.
   */
  int i, j, idxi, idxo, nnus, ncpy, nbytes, sznew, *tdnew, *arr;
  unsigned int d;
  hx_array anew;

  /* determine the number of nonuniform dimensions. */
  for (d = 0, nnus = 0; d < D->nd; d++)
    nnus += (D->dims[d].nus ? 1 : 0);

  /* return successfully if all dimensions are uniform. */
  if (nnus == 0)
    return 1;

  /* check that the datum contains a schedule array. */
  if (D->sched == NULL || D->d_sched < 1 || D->n_sched < 1)
    throw("datum contains no schedule array");

  /* allocate a new size array to build the new array. */
  tdnew = hx_array_index_alloc(D->nd);

  /* check that the size array was allocated. */
  if (!tdnew)
    throw("failed to allocate %u indices", D->nd);

  /* initialize the calculated values. */
  sznew = ncpy = D->dims[0].sz;
  tdnew[0] = D->dims[0].td;

  /* loop over the datum dimensions to build the new size array. */
  for (d = 1; d < D->nd; d++) {
    /* compute the new time-domain point count. */
    tdnew[d] = (D->dims[d].nus ? D->dims[d].tdunif : D->dims[d].td);

    /* compute the number of points to copy per trace. */
    ncpy *= (D->dims[d].cx ? 2 : 1);

    /* compute the array length. */
    sznew *= tdnew[d];
  }

  /* compute the number of bytes to copy per trace. */
  nbytes = ncpy * D->array.n * sizeof(real);

  /* allocate a new array to store the infilled values. */
  if (!hx_array_alloc(&anew, D->array.d, D->array.k, &sznew))
    throw("failed to allocate new (%d, %d)-array", D->array.d, D->array.k);

  /* allocate a new index array for traversal. */
  arr = hx_array_index_alloc(D->nd);

  /* check that allocation was successful. */
  if (!arr)
    throw("failed to allocate %u indices", D->nd);

  /* loop over the indices in the schedule. */
  for (i = 0; i < D->n_sched; i++) {
    /* build the current unpacked index array. */
    for (j = 0, arr[0] = 0; j < D->d_sched; j++)
      arr[j + 1] = D->sched[i * D->d_sched + j];

    /* compute the input array index. */
    idxi = i * ncpy * D->array.n;

    /* compute the output array index. */
    hx_array_index_pack(D->nd, tdnew, arr, &idxo);
    idxo *= anew.n;

    /* copy the current trace into the new array. */
    memcpy(anew.x + idxo, D->array.x + idxi, nbytes);
  }

  /* replaced the datum array with the local infilled array. */
  hx_array_free(&D->array);
  hx_array_copy(&D->array, &anew);
  hx_array_free(&anew);

  /* store the new sizes in the datum dimensions. */
  for (d = 1; d < D->nd; d++) {
    /* store the time-domain size. */
    D->dims[d].sz = D->dims[d].td = tdnew[d];

    /* check if the dimension is complex. */
    if (D->dims[d].cx)
      D->dims[d].sz /= 2;
  }

  /* free the allocated index arrays. */
  free(tdnew);
  free(arr);

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

  /* FIXME: handle cases requiring sign alternation. */

  /* FIXME: handle cases requiring imaginary negation. */

  /* FIXME: handle cases requiring gradient-enhanced arithmetic. */

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

    /* infill nonuniformly sampled indirect dimensions. */
    if (d == 0 && !datum_infill_array(D))
      throw("failed to infill nonuniformly sampled dimensions");
  }

  /* return success. */
  return 1;
}

/* datum_read_sched(): reads a bruker/varian schedule file into the schedule
 * array of a datum structure.
 * @D: pointer to the datum structure to manipulate.
 * @fname: the input filename.
 */
int datum_read_sched (datum *D, const char *fname) {
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
  int d, n, i, *sched;
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
      sched = (int*) realloc(sched, d * n * sizeof(int));

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

  /* correct for group delay *after* the array has been refactored. */
  if (D->type == DATUM_TYPE_BRUKER) {
    /* correct the loaded raw data for bruker group delay uberfail. */
    if (!bruker_fix_grpdelay(&D->array, D->grpdelay))
      throw("failed to correct bruker group delay");

    /* reset the group delay back to zero. */
    D->grpdelay = 0.0;
  }

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

/* datum_resize_array(): resizes each array dimension (if requested) in an
 * NMR datum, but does not change the number of dimensions.
 * @D: pointer to the datum to manipulate.
 * @sz: new size array to use for resizing.
 */
int datum_resize_array (datum *D, int *sz) {
  /* declare a few required variables:
   * @d: dimension loop counter.
   */
  unsigned int d;

  /* check that the array has been allocated. */
  if (!D->array_alloc)
    throw("array is unallocated");

  /* check that each new size is greater than one. */
  for (d = 0; d < D->nd; d++) {
    /* check that the current size is in bounds. */
    if (sz[d] < 2)
      throw("invalid size %d along dimension %u", sz[d], d);
  }

  /* attempt to resize the array content. */
  if (!hx_array_resize(&D->array, D->array.d, D->array.k, sz))
    throw("failed to resize core array");

  /* loop over the acquisition dimensions to store the new sizes. */
  for (d = 0; d < D->nd; d++)
    D->dims[d].sz = sz[d];

  /* return success. */
  return 1;
}

/* datum_slice_array(): slices out a portion of the array from an NMR datum.
 * @D: pointer to the datum to manipulate.
 * @lower: lower bound index array.
 * @upper: upper bound index array.
 */
int datum_slice_array (datum *D, int *lower, int *upper) {
  /* declare a few required variables:
   * @arrnew: destination slice array.
   * @ndnew: new number of dimensions.
   * @ord: dimension reordering array.
   */
  int ndnew, d, dadj, *ord;
  hx_array arrnew;

  /* initialize the local array values, just to be safe. */
  arrnew.d = arrnew.k = 0;
  arrnew.sz = NULL;

  /* slice the datum array into a local array. */
  if (!hx_array_slice(&D->array, &arrnew, lower, upper))
    throw("failed to slice datum array");

  /* replace the datum array with the local array. */
  hx_array_free(&D->array);
  hx_array_copy(&D->array, &arrnew);
  hx_array_free(&arrnew);

  /* allocate an array to specify how to reorder the datum dimensions,
   * if the need arises.
   */
  ord = hx_array_index_alloc(D->nd);

  /* check that allocation succeeded. */
  if (!ord)
    throw("failed to allocate %u indices", D->nd);

  /* store the new array sizes in the datum dimensions. */
  for (d = 0, dadj = 0, ndnew = 0; d < D->nd; d++) {
    /* check if the current dimension has nonzero size. */
    if (D->array.sz[d] > 1) {
      /* store a sortable value in the ordering array. */
      ord[d] = dadj;
      dadj++;

      /* increment the number of nonzero-size dimensions. */
      ndnew++;
    }
    else {
      /* store an unsortable value in the ordering array. */
      ord[d] = (int) D->nd;
    }

    /* store the new array dimension size. */
    D->dims[d].sz = D->array.sz[d];
  }

  /* compact zero-size array dimensions out of the array. */
  if (!hx_array_compact(&D->array))
    throw("failed to compact core array");

  /* check if the dimension count changed. */
  if (ndnew != D->nd) {
    /* reorder the datum dimensions. */
    if (!datum_reorder_dims(D, ord))
      throw("failed to reorder datum dimensions");

    /* reallocate the dimension array. */
    D->dims = (datum_dim*) realloc(D->dims, ndnew * sizeof(datum_dim));

    /* check that allocated succeeded. */
    if (D->dims == NULL)
      throw("failed to reallocate dimension array");

    /* sort the ordering array into a valid index list. */
    hx_array_index_sort(D->nd, ord);

    /* reorder the core array basis elements. */
    if (!hx_array_reorder_bases(&D->array, ord))
      throw("failed to reorder core array basis elements");

    /* reduce the dimensionality of the core array. */
    if (!hx_array_resize(&D->array, ndnew, D->array.k, D->array.sz))
      throw("failed to resize core array");

    /* set the new dimension count. */
    D->nd = ndnew;
  }

  /* free the allocated index array. */
  free(ord);

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

