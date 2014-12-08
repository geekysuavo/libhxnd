
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

/* include the n-dimensional math header. */
#include <hxnd/hx.h>

/* include the byte-level data header. */
#include <hxnd/bytes.h>

/* define the size of the data buffer for raw file reading, in bytes.
 */
#define HX_ARRAY_FREAD_SZ_BUF  33554432

/* hx_array_fread_raw(): read raw bytes from a file out of a given format
 * and into a linear hypercomplex array.
 * @fh: input file handle.
 * @x: output array structure pointer.
 * @endian: raw word byte ordering.
 * @wordsz: number of bytes per word.
 * @isflt: whether the words are floats.
 * @offhead: file header offset in bytes.
 * @offblk: file block offset in bytes.
 * @nblks: number of data blocks.
 * @nwords: number of data words per block.
 * @nalign: number of bytes to align blocks to.
 */
int hx_array_fread_raw (FILE *fh, hx_array *x, enum byteorder endian,
                        unsigned int wordsz, unsigned int isflt,
                        unsigned int offhead, unsigned int offblk,
                        unsigned int nblks, unsigned int nwords,
                        unsigned int nalign) {
  /* declare variables required for buffered reading:
   * @pos: current file position, in bytes.
   * @i: main block loop counter.
   * @n: number of blocks to read.
   * @nbytes: number of bytes per block.
   * @nbuf: number of buffer bytes.
   * @nask: number to requested bytes.
   * @nread: number of read bytes.
   * @nrem: number of bytes remaining.
   * @ialign: number of alignment offsets.
   * @buf: array of buffer bytes.
   */
  unsigned int pos, i, n, nbytes, nbuf, nask, nread, nrem, ialign;
  uint8_t *buf;

  /* declare variables required for coefficient allocation and access:
   * @k: raw buffer byte index for conversion into array words.
   * @xi: current array coefficient pointer.
   * @len: total number of words.
   */
  unsigned int k;
  real *xi;
  int len;

  /* compute the number of bytes per block, and the total byte count. */
  nbytes = nwords * wordsz;
  nrem = nbytes * nblks;

  /* compute the total number of data words. */
  len = nblks * nwords;

  /* initialize the block count and size. */
  nbuf = nbytes;
  n = nblks;

  /* determine whether block alignment is specified. */
  if (nalign) {
    /* check for the special case in which the reads may be lumped:
     * 1. the blocks fall on the alignment boundaries.
     * 2. no block header exists.
     * 3. no file header exists.
     */
    if (offhead == 0 && offblk == 0 && nbytes % nalign == 0) {
      /* in this special case, we can increase the buffer size. */
      nbuf *= n;
      n = 1;
    }
  }

  /* determine whether logical blocking must be done. */
  if (n == 1 && nbuf > HX_ARRAY_FREAD_SZ_BUF) {
    /* logical blocking is needed. determine the number of logical blocks:
     *  - if the buffer is evenly divided by the max size, use that ratio.
     *  - if not, add one to the ratio, where the last read will be partial.
     */
    n = (nbuf % HX_ARRAY_FREAD_SZ_BUF ?
         nbuf / HX_ARRAY_FREAD_SZ_BUF + 1 :
         nbuf / HX_ARRAY_FREAD_SZ_BUF);

    /* use the maximum allowed read buffer size. */
    nbuf = HX_ARRAY_FREAD_SZ_BUF;
  }

  /* allocate memory for the buffer bytes. */
  buf = (uint8_t*) calloc(nbuf, sizeof(uint8_t));

  /* check that memory was allocated. */
  if (!buf)
    throw("failed to allocate %u buffer bytes", nbytes);

  /* allocate the destination array structure. */
  if (!hx_array_alloc(x, 0, 1, &len))
    throw("failed to allocate %d array coefficients", len);

  /* move past the file header. */
  if (offhead && fseek(fh, offhead, SEEK_SET))
    throw("failed to seek %u bytes past file header", offhead);

  /* initialize the file position and the alignment offset. */
  pos = offhead;
  ialign = 0;

  /* loop over the data blocks. */
  for (i = 0, xi = x->x; i < n; i++) {
    /* check if alignment is required. */
    if (nalign) {
      /* align to the next boundary. */
      ialign = pos / nalign;
      while (pos % nalign)
        pos = ++ialign * nalign;

      /* move to the new aligned location. */
      if (fseek(fh, pos, SEEK_SET))
        throw("failed to seek to %u-byte alignment boundary", nalign);
    }

    /* move past the block header. */
    pos += offblk;
    if (offblk && fseek(fh, offblk, SEEK_CUR))
      throw("failed to seek %u bytes past block header", offblk);

    /* determine how many bytes to ask for. */
    nask = (nrem > nbuf ? nbuf : nrem);

    /* read the current data block. */
    nread = fread(buf, sizeof(uint8_t), nask, fh);

    /* check that the read succeeded. */
    if (nread != nask)
      throw("failed to read data block #%u", i);

    /* subtract the read bytes from the bytes remaining. */
    nrem -= nread;

    /* check if byte swaps are required. */
    if (!bytes_native(endian) && wordsz > 1)
      bytes_swap(buf, nbuf / wordsz, wordsz);

    /* copy the read words into the final array. */
    for (k = 0; k < nread; k += wordsz)
      *(xi++) = bytes_unpack(buf + k, wordsz, isflt);
  }

  /* free the allocated buffer memory. */
  free(buf);

  /* return success. */
  return 1;
}

/* hx_array_fwrite_raw(): write a hypercomplex array into a raw output format
 * that differs from the internal storage format of the hx_array.
 * @fh: output file handle.
 * @x: source array structure pointer.
 * @endian: raw word byte ordering.
 * @wordsz: number of bytes per word.
 * @isflt: whether the words are floats.
 */
int hx_array_fwrite_raw (FILE *fh, hx_array *x, enum byteorder endian,
                         unsigned int wordsz, unsigned int isflt) {
  /* declare a few required variables:
   * @bytes: raw array of bytes, size of a word.
   * @i: general-purpose loop counter.
   */
  uint8_t *bytes;
  int i;

  /* check if the requested parameters match the internal storage. */
  if (bytes_native(endian) && wordsz == sizeof(real) && isflt) {
    /* use a much faster single write. */
    if (fwrite(x->x, sizeof(real), x->len, fh) != x->len)
      throw("failed to write %d %cE %u-byte %ss", x->len,
            endian == BYTES_ENDIAN_LITTLE ? 'L' :
            endian == BYTES_ENDIAN_BIG ? 'B' : 'U',
            wordsz, isflt ? "float" : "integer");

    /* return successfully. */
    return 1;
  }

  /* allocate the raw byte array. */
  bytes = (uint8_t*) malloc(wordsz * sizeof(uint8_t));

  /* check that allocation was successful. */
  if (!bytes)
    throw("failed to allocate a %u-byte storage array", wordsz);

  /* loop over the array coefficients. */
  for (i = 0; i < x->len; i++) {
    /* pack the coefficient into the byte array. */
    if (!bytes_pack(x->x[i], bytes, wordsz, isflt))
      throw("failed to pack coefficient %d", i);

    /* write the raw bytes out to the file. */
    if (fwrite(bytes, sizeof(uint8_t), wordsz, fh) != 1)
      throw("failed to write %cE %u-byte %s #%d",
            endian == BYTES_ENDIAN_LITTLE ? 'L' :
            endian == BYTES_ENDIAN_BIG ? 'B' : 'U',
            wordsz, isflt ? "float" : "integer", i);
  }

  /* free the raw byte array. */
  free(bytes);

  /* return success. */
  return 1;
}

