
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

/* include the pipe header. */
#include <hxnd/nmr-pipe.h>

/* pipe_check_magic(): checks the first "magic" bytes of a file and
 * returns whether they match the pipe-format file magic number.
 * @fname: the input filename.
 */
int pipe_check_magic (const char *fname) {
  /* declare a few required variables:
   * @fh: input file handle.
   * @wd: first words of input file.
   */
  float wd[3];
  FILE *fh;

  /* open the input file. */
  fh = fopen(fname, "rb");

  /* check that the file was opened. */
  if (!fh)
    throw("failed to open '%s'", fname);

  /* read the first three words. */
  if (!fread(&wd, sizeof(float), 3, fh))
    throw("failed to read magic numbers");

  /* close the input file. */
  fclose(fh);

  /* check the magic word, without swapping. */
  if (wd[2] == (float) PIPE_MAGIC)
    return 1;

  /* swap the word and check again. */
  bytes_swap_u32((uint32_t*) &wd[2]);
  if (wd[2] == (float) PIPE_MAGIC)
    return 1;

  /* no match. */
  return 0;
}

/* pipe_read_header(): read the contents of a header from a pipe-format file.
 * @fname: the input filename.
 * @endianness: byte order result pointer.
 * @hdr: header result pointer.
 */
int pipe_read_header (const char *fname, enum byteorder *endianness,
                      struct pipe_header *hdr) {
  /* declare a few required variables:
   */
  unsigned int n_bytes;
  uint8_t *bytes;

  /* read in the file header bytes. */
  n_bytes = sizeof(struct pipe_header);
  bytes = bytes_read_block(fname, 0, n_bytes);

  /* check that the bytes were read successfully. */
  if (!bytes)
    throw("failed to read header from '%s'", fname);

  /* copy the header bytes onto the header structure. */
  memcpy(hdr, bytes, n_bytes);

  /* check if the floating point order value is correct. if not, the bytes
   * likely need to be swapped.
   */
  if (hdr->order != 0.0 && hdr->order != (float) PIPE_MAGIC) {
    /* swap the bytes of each word. */
    bytes_swap_general(bytes, n_bytes, sizeof(float));

    /* re-copy the bytes onto the structure. */
    memcpy(hdr, bytes, n_bytes);

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

  /* free the read bytes. */
  free(bytes);

  /* return success. */
  return 1;
}

/* pipe_print_header(): prints the contents of a pipe file header to the
 * standard output stream.
 * @hdr: pointer to the header to print.
 */
void pipe_print_header (struct pipe_header *hdr) {
  /* print the global values. */
  printf("magic, format, order = %.1f, %.0f, %.3f\n",
         hdr->magic, hdr->format, hdr->order);

  /* print the dimension count. */
  printf("dims = %.0f (%stransposed)\n", hdr->ndims,
         hdr->trans ? "" : "not ");

  /* print the dimension ordering. */
  printf("dim order = %.0f, %.0f, %.0f, %.0f\n",
         hdr->dimorder[0], hdr->dimorder[1],
         hdr->dimorder[2], hdr->dimorder[3]);

  /* print the dimension-specific values header. */
  printf("%-15s%15s%15s%15s%15s\n", "", "f1", "f2", "f3", "f4");

  /* print the dimension labels. */
  printf("%-15s%15s%15s%15s%15s\n", "",
         hdr->label_f1, hdr->label_f2,
         hdr->label_f3, hdr->label_f4);

  /* print the apodization point counts. */
  printf("%-15s%15.0f%15.0f%15.0f%15.0f\n", "apod",
         hdr->apod_f1, hdr->apod_f2,
         hdr->apod_f3, hdr->apod_f4);

  /* print the spectral widths. */
  printf("%-15s%15.3f%15.3f%15.3f%15.3f\n", "sw",
         hdr->sw_f1, hdr->sw_f2,
         hdr->sw_f3, hdr->sw_f4);

  /* print the observe frequencies. */
  printf("%-15s%15.3f%15.3f%15.3f%15.3f\n", "obs",
         hdr->obs_f1, hdr->obs_f2,
         hdr->obs_f3, hdr->obs_f4);

  /* print the origin frequencies. */
  printf("%-15s%15.3f%15.3f%15.3f%15.3f\n", "orig",
         hdr->orig_f1, hdr->orig_f2,
         hdr->orig_f3, hdr->orig_f4);

  /* print the dimension sizes. */
  printf("%-15s%15.0f%15.0f%15.0f%15.0f\n", "size",
         hdr->sz, hdr->specnum,
         hdr->size_f3, hdr->size_f4);

  /* print the time-domain point counts. */
  printf("%-15s%15.0f%15.0f%15.0f%15.0f\n", "tdsz",
         hdr->tdsz_f1, hdr->tdsz_f2,
         hdr->tdsz_f3, hdr->tdsz_f4);

  /* print the zero-fill size. */
  printf("%-15s%15.0f%15.0f%15.0f%15.0f\n", "zf",
         hdr->zf_f1, hdr->zf_f2,
         hdr->zf_f3, hdr->zf_f4);

  /* print the fft point counts. */
  printf("%-15s%15.0f%15.0f%15.0f%15.0f\n", "ftsz",
         hdr->ftsz_f1, hdr->ftsz_f2,
         hdr->ftsz_f3, hdr->ftsz_f4);

  /* print the extracted point counts. */
  printf("%-15s%15.0f%15.0f%15.0f%15.0f\n", "x1",
         hdr->x1_f1, hdr->x1_f2,
         hdr->x1_f3, hdr->x1_f4);

  /* print the extracted point counts. */
  printf("%-15s%15.0f%15.0f%15.0f%15.0f\n", "xn",
         hdr->xn_f1, hdr->xn_f2,
         hdr->xn_f3, hdr->xn_f4);

  /* print the quadrature flags. */
  printf("%-15s%15.0f%15.0f%15.0f%15.0f\n", "quad",
         hdr->quad_f1, hdr->quad_f2,
         hdr->quad_f3, hdr->quad_f4);

  /* print the fft flags. */
  printf("%-15s%15.0f%15.0f%15.0f%15.0f\n", "ft",
         hdr->ftflag_f1, hdr->ftflag_f2,
         hdr->ftflag_f3, hdr->ftflag_f4);

  /* print the acquisition sign flags. */
  printf("%-15s%15.0f%15.0f%15.0f%15.0f\n", "aqsgn",
         hdr->aqsgn_f1, hdr->aqsgn_f2,
         hdr->aqsgn_f3, hdr->aqsgn_f4);
}

/* pipe_read(): reads a pipe-format data file into a real linear array.
 * @fname: the input data filename.
 * @n: the number of expected data bytes.
 * @x: the output array.
 */
int pipe_read (const char *fname, unsigned int n, hx_array *x) {
  /* declare a few required variables. */
  enum byteorder endianness = BYTES_ENDIAN_AUTO;
  unsigned int n_actual, offset;
  struct pipe_header hdr;
  uint8_t *bytes;

  /* read the header information from the data file. */
  if (!pipe_read_header(fname, &endianness, &hdr))
    throw("failed to read header of '%s'", fname);

  /* compute the byte offset from which to begin reading point data. */
  offset = sizeof(struct pipe_header);

  /* compute the actual number of bytes in the file. */
  n_actual = bytes_size(fname) - offset;

  /* check that the byte counts match. */
  if (n != n_actual)
    throw("data size mismatch (expected %u, read %u)", n, n_actual);

  /* read the data bytes in from the file. */
  bytes = bytes_read_block(fname, offset, n);

  /* check that the bytes were read successfully. */
  if (!bytes)
    throw("failed to read %u bytes from '%s'", n, fname);

  /* build a real linear array from the byte data. */
  if (!bytes_toarray(bytes, n, endianness, sizeof(float), 1, x))
    throw("failed to convert bytes to array");

  /* free the read byte data. */
  free(bytes);

  /* return success. */
  return 1;
}

/* pipe_interlace(): pipe-format files do not store their individual traces
 * as interlaced R,I,R,I..., but instead interlaces entire real traces and
 * imaginary traces. this fixes pipe's mistake so hx_array_deinterlace()
 * can operate on pipe data just like any other data type.
 * @x: pointer to the array to interlace.
 * @n: number of *complex* points per trace.
 */
int pipe_interlace (hx_array *x, unsigned int n) {
  /* declare a few required variables:
   * @i: trace loop counter.
   * @j: point loop counter.
   * @ntraces: number of complex traces.
   * @kre: real trace start index.
   * @kim: imag trace start index.
   * @xtmp: temporary duplicate coefficient array.
   */
  unsigned int i, j, ntraces, kre, kim;
  real *xtmp;

  /* check that the trace size evenly divides the array elements. */
  if (x->sz[0] % (2 * n))
    throw("trace size %u does not evenly divide array (%d)", n, x->sz[0]);

  /* compute the number of traces in the array. */
  ntraces = x->sz[0] / (2 * n);

  /* allocate a temporary duplicate array of the array elements. */
  xtmp = (real*) malloc(x->len * sizeof(real));
  if (!xtmp)
    throw("failed to allocate temporary array of %d reals", x->len);

  /* copy the data over into the duplicate array. */
  memcpy(xtmp, x->x, x->len * sizeof(real));

  /* loop over each trace in the linear real array. */
  for (i = 0; i < ntraces; i++) {
    /* compute the real and imaginary trace starting indices. */
    kre = (2 * i) * n;
    kim = (2 * i + 1) * n;

    /* loop over the points of the trace. */
    for (j = 0; j < n; j++) {
      /* shuffle the points around. */
      x->x[i * n + 2 * j] = xtmp[kre + j];
      x->x[i * n + 2 * j + 1] = xtmp[kim + j];
    }
  }

  /* free the duplicate array. */
  free(xtmp);

  /* return success. */
  return 1;
}

/* pipe_fill_datum(): intelligently parses pipe parameters into an
 * NMR datum structure.
 * @fname: the input filename.
 * @D: pointer to the datum struct to fill.
 */
int pipe_fill_datum (const char *fname, datum *D) {
  /* declare variables required to determine byte ordering:
   * @endianness: the byte ordering of the data file.
   * @hdr: the pipe file header structure.
   */
  enum byteorder endianness = BYTES_ENDIAN_AUTO;
  struct pipe_header hdr;

  /* declare variables required to traverse dimensions:
   * @d: dimension loop counter.
   * @ord: dimension index array.
   */
  unsigned int d;
  int ord[PIPE_MAXDIM];

  /* read the header information from the data file. */
  if (!pipe_read_header(fname, &endianness, &hdr))
    throw("failed to read header of '%s'", fname);

  /* initially set the number of dimensions to the maximum allowed, because
   * pipe arranges its dimension information in a really screwy way.
   */
  D->nd = PIPE_MAXDIM;

  /* check the dimensionality. */
  if ((int) hdr.ndims < 1 ||
      (int) hdr.ndims > PIPE_MAXDIM)
    throw("invalid dimensionality %.0f", hdr.ndims);

  /* allocate the dimension parameter array. */
  D->dims = (datum_dim*) calloc(PIPE_MAXDIM, sizeof(datum_dim));

  /* check that the dimension parameter array was allocated. */
  if (D->dims == NULL)
    throw("failed to allocate four datum dimensions");

  /* store the dimension ordering array values. */
  for (d = 0; d < PIPE_MAXDIM; d++)
    ord[d] = hdr.dimorder[d] - 1;

  /* store the quadrature flags. */
  D->dims[ord[0]].cx = ((int) hdr.quad_f1 != 1);
  D->dims[ord[1]].cx = ((int) hdr.quad_f2 != 1);
  D->dims[ord[2]].cx = ((int) hdr.quad_f3 != 1);
  D->dims[ord[3]].cx = ((int) hdr.quad_f4 != 1);

  /* store the fourier-transform flags. */
  D->dims[ord[0]].ft = (unsigned int) hdr.ftflag_f1;
  D->dims[ord[1]].ft = (unsigned int) hdr.ftflag_f2;
  D->dims[ord[2]].ft = (unsigned int) hdr.ftflag_f3;
  D->dims[ord[3]].ft = (unsigned int) hdr.ftflag_f4;

  /* store the nucleus strings. */
  strncpy(D->dims[ord[0]].nuc, hdr.label_f1, PIPE_HDRSTR_SZ_LABEL);
  strncpy(D->dims[ord[1]].nuc, hdr.label_f2, PIPE_HDRSTR_SZ_LABEL);
  strncpy(D->dims[ord[2]].nuc, hdr.label_f3, PIPE_HDRSTR_SZ_LABEL);
  strncpy(D->dims[ord[3]].nuc, hdr.label_f4, PIPE_HDRSTR_SZ_LABEL);

  /* null-terminate the nucleus string. */
  for (d = 0; d < PIPE_MAXDIM; d++)
    D->dims[d].nuc[7] = '\0';

  /* store the size parameters. */
  D->dims[ord[0]].td = D->dims[ord[0]].tdunif = hdr.tdsz_f1;
  D->dims[ord[1]].td = D->dims[ord[1]].tdunif = hdr.tdsz_f2;
  D->dims[ord[2]].td = D->dims[ord[2]].tdunif = hdr.tdsz_f3;
  D->dims[ord[3]].td = D->dims[ord[3]].tdunif = hdr.tdsz_f4;

  /* store the f1 current size parameter. */
  D->dims[ord[0]].sz = hdr.x1_f1 && hdr.xn_f1 ? hdr.xn_f1 - hdr.x1_f1 + 1 :
                       D->dims[ord[0]].ft ? hdr.ftsz_f1 : hdr.apod_f1;

  /* store the f2 current size parameter. */
  D->dims[ord[1]].sz = hdr.x1_f2 && hdr.xn_f2 ? hdr.xn_f2 - hdr.x1_f2 + 1 :
                       D->dims[ord[1]].ft ? hdr.ftsz_f2 : hdr.apod_f2;

  /* store the f3 current size parameter. */
  D->dims[ord[2]].sz = hdr.x1_f3 && hdr.xn_f3 ? hdr.xn_f3 - hdr.x1_f3 + 1 :
                       D->dims[ord[2]].ft ? hdr.ftsz_f3 : hdr.apod_f3;

  /* store the f4 current size parameter. */
  D->dims[ord[3]].sz = hdr.x1_f4 && hdr.xn_f4 ? hdr.xn_f4 - hdr.x1_f4 + 1 :
                       D->dims[ord[3]].ft ? hdr.ftsz_f4 : hdr.apod_f4;

  /* store the spectral width parameters. */
  D->dims[ord[0]].width = hdr.sw_f1;
  D->dims[ord[1]].width = hdr.sw_f2;
  D->dims[ord[2]].width = hdr.sw_f3;
  D->dims[ord[3]].width = hdr.sw_f4;

  /* store the carrier frequency parameters. */
  D->dims[ord[0]].carrier = hdr.obs_f1;
  D->dims[ord[1]].carrier = hdr.obs_f2;
  D->dims[ord[2]].carrier = hdr.obs_f3;
  D->dims[ord[3]].carrier = hdr.obs_f4;

  /* store the spectral offset parameters. */
  D->dims[ord[0]].offset = hdr.car_f1 * hdr.obs_f1;
  D->dims[ord[1]].offset = hdr.car_f2 * hdr.obs_f2;
  D->dims[ord[2]].offset = hdr.car_f3 * hdr.obs_f3;
  D->dims[ord[3]].offset = hdr.car_f4 * hdr.obs_f4;

  /* set the true dimension count and reallocate the dimension array. */
  D->nd = (unsigned int) hdr.ndims;
  D->dims = (datum_dim*) realloc(D->dims, D->nd * sizeof(datum_dim));

  /* check that the array was reallocated successfully. */
  if (D->dims == NULL)
    throw("failed to resize dimension array to %u elements", D->nd);

  /* store the filename string. */
  D->fname = (char*) malloc((strlen(fname) + 1) * sizeof(char));
  if (D->fname)
    strcpy(D->fname, fname);

  /* store the datum type. */
  D->type = DATUM_TYPE_PIPE;
  D->endian = endianness;

  /* return success. */
  return 1;
}

