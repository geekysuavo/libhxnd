
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

/* include the hx native format header. */
#include <hxnd/nmr-hxnd.h>

/* define the number of (u64) members in the header and dimension sections of
 * binary datum files.
 */
#define HXND_SZ_HDR   7
#define HXND_SZ_DIM  10

/* define bit field positions to store status flags in binary datum files.
 */
#define HXND_S_COMPLEX  0x0000000000000001
#define HXND_S_NUS      0x0000000000000002
#define HXND_S_FFT      0x0000000000000004
#define HXND_S_ALT      0x0000000000000008
#define HXND_S_NEG      0x0000000000000010
#define HXND_S_GENH     0x0000000000000020

/* hxnd_guess(): check whether a file contains hx-native nmr data.
 * @fname: the input filename.
 */
int hxnd_guess (const char *fname) {
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
  if (wd == HXND_MAGIC)
    return 1;

  /* swap the word and check again. */
  bytes_swap_u64(&wd);
  if (wd == HXND_MAGIC)
    return 1;

  /* no match. */
  return 0;
}

/* hxnd_decode(): read datum parameters from a file.
 * @D: pointer to the destination datum structure.
 * @fname: the input filename.
 */
int hxnd_decode (datum *D, const char *fname) {
  /* declare a few required variables:
   * @i: buffer index.
   * @d: dimension loop counter.
   * @n_buf: current buffer word size.
   * @buf: output header/dimension buffer.
   * @status: status word value.
   * @fh: input file handle.
   */
  unsigned int i, d, n_buf, n_sched, swapping;
  uint64_t *buf, status;
  FILE *fh;

  /* open the input file. */
  if (fname)
    fh = fopen(fname, "rb");
  else
    fh = stdin;

  /* check that the file was opened. */
  if (!fh)
    throw("failed to open '%s'", fname);

  /* allocate the header buffer. */
  n_buf = HXND_SZ_HDR;
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
  if (buf[0] != HXND_MAGIC) {
    /* no match. swap the bytes of each word. */
    bytes_swap((uint8_t*) buf, n_buf, sizeof(uint64_t));

    /* now check the magic word. */
    if (buf[0] != HXND_MAGIC)
      throw("invalid magic number 0x%016lx", buf[0]);

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
  D->epoch = (time_t) buf[i++];
  D->nd = (unsigned int) buf[i++];
  D->d_sched = (int) buf[i++];
  D->n_sched = (int) buf[i++];

  /* check if any schedule was stored in the file. */
  n_sched = D->n_sched * D->d_sched;
  if (n_sched) {
    /* allocate the schedule array. */
    D->sched = hx_array_index_alloc(n_sched);

    /* check that allocation succeeded. */
    if (D->sched == NULL)
      throw("failed to allocate %u schedule indices", n_sched);

    /* reallocate the buffer to hold schedule values. */
    n_buf = n_sched;
    buf = (uint64_t*) realloc(buf, n_buf * sizeof(uint64_t));

    /* check that allocation succeeded. */
    if (!buf)
      throw("failed to allocate schedule buffer");

    /* read the schedule data. */
    if (fread(buf, sizeof(uint64_t), n_buf, fh) != n_buf)
      throw("failed to read schedule");

    /* swap the bytes, if required. */
    if (swapping)
      bytes_swap((uint8_t*) buf, n_buf, sizeof(uint64_t));

    /* unpack the schedule values. */
    for (i = 0; i < n_sched; i++)
      D->sched[i] = (int) buf[i];
  }
  else {
    /* null the schedule. */
    D->sched = NULL;
  }

  /* allocate the dimension array. */
  D->dims = (datum_dim*) calloc(D->nd, sizeof(datum_dim));

  /* check that allocation succeeded. */
  if (D->dims == NULL)
    throw("failed to allocate %u dimensions", D->nd);

  /* reallocate the buffer to hold dimension information. */
  n_buf = HXND_SZ_DIM;
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

    /* unpack the array dimension indices. */
    i = 0;
    D->dims[d].d = *((int*) (buf + (i++)));
    D->dims[d].k = *((int*) (buf + (i++)));

    /* unpack the dimension size parameters. */
    D->dims[d].sz = (unsigned int) buf[i++];
    D->dims[d].td = (unsigned int) buf[i++];
    D->dims[d].tdunif = (unsigned int) buf[i++];

    /* unpack the status word. */
    status = buf[i++];
    D->dims[d].cx = (status &   HXND_S_COMPLEX ? 1 : 0);
    D->dims[d].nus = (status &  HXND_S_NUS ? 1 : 0);
    D->dims[d].ft = (status &   HXND_S_FFT ? 1 : 0);
    D->dims[d].alt = (status &  HXND_S_ALT ? 1 : 0);
    D->dims[d].neg = (status &  HXND_S_NEG ? 1 : 0);
    D->dims[d].genh = (status & HXND_S_GENH ? 1 : 0);

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
  if (fname == NULL) {
    /* read the core array content from the end of the file stream. */
    if (!hx_array_fread(&D->array, fh))
      throw("failed to read core array");

    /* read succeeded. set the array allocation status flag. */
    D->array_alloc = 1;
  }

  /* close the input file and store the filename in the datum. */
  if (fname) {
    /* close the file. */
    fclose(fh);

    /* store the filename in the datum. */
    D->fname = (char*) malloc((strlen(fname) + 1) * sizeof(char));
    if (D->fname)
      strcpy(D->fname, fname);
  }

  /* return success. */
  return 1;
}

/* hxnd_encode(): write a datum structure to a file.
 * @D: pointer to the source datum structure.
 * @fname: the output filename.
 */
int hxnd_encode (datum *D, const char *fname) {
  /* declare a few required variables:
   * @i: buffer index.
   * @d: dimension loop counter.
   * @n_buf: current buffer word size.
   * @buf: output header/dimension buffer.
   * @status: status word value.
   * @fh: output file handle.
   */
  unsigned int i, j, d, n_buf, n_sched;
  uint64_t *buf, status;
  FILE *fh;

  /* open the output file. */
  if (fname)
    fh = fopen(fname, "wb");
  else
    fh = stdout;

  /* check that the file was opened. */
  if (!fh)
    throw("failed to open '%s'", fname);

  /* allocate the header buffer. */
  n_sched = D->n_sched * D->d_sched;
  n_buf = HXND_SZ_HDR + n_sched;
  buf = (uint64_t*) calloc(n_buf, sizeof(uint64_t));

  /* check that allocation succeeded. */
  if (!buf)
    throw("failed to allocate header buffer");

  /* build the file header. */
  i = 0;
  buf[i++] = (uint64_t) HXND_MAGIC;
  buf[i++] = (uint64_t) D->endian;
  buf[i++] = (uint64_t) DATUM_TYPE_HXND;
  buf[i++] = (uint64_t) D->epoch;
  buf[i++] = (uint64_t) D->nd;
  buf[i++] = (uint64_t) D->d_sched;
  buf[i++] = (uint64_t) D->n_sched;

  /* add the schedule values into the header. */
  for (j = 0; j < n_sched; j++)
    buf[i++] = (uint64_t) D->sched[j];

  /* write the file header. */
  if (fwrite(buf, sizeof(uint64_t), n_buf, fh) != n_buf)
    throw("failed to write header");

  /* reallocate the buffer to hold dimension information. */
  n_buf = HXND_SZ_DIM;
  buf = (uint64_t*) realloc(buf, n_buf * sizeof(uint64_t));

  /* check that reallocation succeeded. */
  if (!buf)
    throw("failed to allocate dimension buffer");

  /* loop over the dimensions in the datum structure. */
  for (d = 0; d < D->nd; d++) {
    /* zero the buffer memory. */
    memset(buf, 0, n_buf * sizeof(uint64_t));

    /* store the array dimension indices. */
    i = 0;
    buf[i++] = (uint64_t) *((unsigned int*) &D->dims[d].d);
    buf[i++] = (uint64_t) *((unsigned int*) &D->dims[d].k);

    /* store the size parameters. */
    buf[i++] = (uint64_t) D->dims[d].sz;
    buf[i++] = (uint64_t) D->dims[d].td;
    buf[i++] = (uint64_t) D->dims[d].tdunif;

    /* build the status word. */
    status = 0;
    status |= (D->dims[d].cx ?   HXND_S_COMPLEX : 0);
    status |= (D->dims[d].nus ?  HXND_S_NUS : 0);
    status |= (D->dims[d].ft ?   HXND_S_FFT : 0);
    status |= (D->dims[d].alt ?  HXND_S_ALT : 0);
    status |= (D->dims[d].neg ?  HXND_S_NEG : 0);
    status |= (D->dims[d].genh ? HXND_S_GENH : 0);

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

  /* close the output file. */
  if (fname)
    fclose(fh);

  /* return success. */
  return 1;
}

/* hxnd_array(): read hx-native array data into a datum structure.
 * @D: pointer to the destination datum structure.
 */
int hxnd_array (datum *D) {
  /* declare a few required variables:
   * @offset: byte offset designating the start of point data.
   * @fh: input file handle.
   */
  unsigned int offset;
  FILE *fh;

  /* open the input file. */
  if (D->fname)
    fh = fopen(D->fname, "rb");
  else
    throw("invalid input filename");

  /* compute the offset. */
  offset = HXND_SZ_HDR;
  offset += D->d_sched * D->n_sched;
  offset += D->nd * HXND_SZ_DIM;
  offset *= sizeof(uint64_t);

  /* seek past the header bytes. */
  if (fseek(fh, offset, SEEK_SET))
    throw("failed to seek to array in '%s'", D->fname);

  /* read the array data from the hx-format file. */
  if (!hx_array_fread(&D->array, fh))
    throw("failed to read array from '%s'", D->fname);

  /* close the input file. */
  fclose(fh);

  /* return success. */
  return 1;
}

