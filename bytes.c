
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

/* include the bytes header. */
#include "bytes.h"

/* bytes_swap_u16(): swaps the bytes of a two-byte word.
 * @x: pointer to the word to swap in-place.
 */
void bytes_swap_u16 (uint16_t *x) {
  /* extract the bytes of the word. */
  uint8_t b0 = BYTES_BYTE0(*x);
  uint8_t b1 = BYTES_BYTE1(*x);

  /* rebuild the word with swapped bytes. */
  *x = BYTES_U16(b1, b0);
}

/* bytes_swap_u32(): swaps the bytes of a four-byte word.
 * @x: pointer to the word to swap in-place.
 */
void bytes_swap_u32 (uint32_t *x) {
  /* extract the bytes of the word. */
  uint8_t b0 = BYTES_BYTE0(*x);
  uint8_t b1 = BYTES_BYTE1(*x);
  uint8_t b2 = BYTES_BYTE2(*x);
  uint8_t b3 = BYTES_BYTE3(*x);

  /* rebuild the word with swapped bytes. */
  *x = BYTES_U32(b3, b2, b1, b0);
}

/* bytes_swap_u64(): swaps the bytes of an eight-byte word.
 * @x: pointer to the word to swap in-place.
 */
void bytes_swap_u64 (uint64_t *x) {
  /* extract the first four bytes of the word. */
  uint8_t b0 = BYTES_BYTE0(*x);
  uint8_t b1 = BYTES_BYTE1(*x);
  uint8_t b2 = BYTES_BYTE2(*x);
  uint8_t b3 = BYTES_BYTE3(*x);

  /* extract the next four bytes of the word. */
  uint32_t y = BYTES_RSHIFT4(*x);
  uint8_t b4 = BYTES_BYTE0(y);
  uint8_t b5 = BYTES_BYTE1(y);
  uint8_t b6 = BYTES_BYTE2(y);
  uint8_t b7 = BYTES_BYTE3(y);

  /* rebuild the word with swapped bytes. */
  *x = BYTES_U64(b7, b6, b5, b4, b3, b2, b1, b0);
}

/* bytes_size(): read the number of bytes in a specified file.
 * @fname: the input filename.
 */
unsigned int bytes_size (const char *fname) {
  /* declare a few required variables. */
  unsigned int n;
  FILE *fh;

  /* open the input file. */
  fh = fopen(fname, "rb");

  /* check that the file was opened successfully. */
  if (!fh)
    return 0;

  /* move to the end of the file. */
  if (fseek(fh, 0, SEEK_END))
    return 0;

  /* read the file size. */
  n = (unsigned int) ftell(fh);

  /* move back to the beginning of the file. */
  if (fseek(fh, 0, SEEK_SET))
    return 0;

  /* close the input file and return the byte count. */
  fclose(fh);
  return n;
}

/* bytes_read_block(): read a specified number of bytes from a binary file,
 * starting at a given offset.
 * @fname: the input filename.
 * @offset: byte offset from the beginning of the file.
 * @n: number of bytes to read.
 */
uint8_t *bytes_read_block (const char *fname,
                           unsigned int offset,
                           unsigned int n) {
  /* declare a few required variables. */
  uint8_t *bytes;
  FILE *fh;

  /* open the input file. */
  fh = fopen(fname, "rb");

  /* check that the file was opened successfully. */
  if (!fh)
    return NULL;

  /* allocate memory for the read bytes. */
  bytes = (uint8_t*) malloc(n * sizeof(uint8_t));

  /* check that the memory was allocated. */
  if (!bytes)
    return NULL;

  /* move back to the beginning of the file. */
  if (fseek(fh, offset, SEEK_SET)) {
    /* free the byte array and return nothing. */
    free(bytes);
    return NULL;
  }

  /* read the bytes in from the input file. */
  if (!fread(bytes, sizeof(uint8_t), n, fh)) {
    /* free the byte array and return nothing. */
    free(bytes);
    return NULL;
  }

  /* close the input file and return the read data. */
  fclose(fh);
  return bytes;
}

/* bytes_read_bruker(): reads fid byte data from a bruker fid/ser file.
 * @fname: the input filename.
 * @nblk: the number of fids in the file.
 * @szblk: the number of bytes per fid.
 * @n: pointer to the output byte count.
 */
uint8_t *bytes_read_bruker (const char *fname,
                            unsigned int nblk,
                            unsigned int szblk,
                            unsigned int *n) {
  /* declare a few required variables:
   * @ndata: number of bytes of real data in the file.
   * @noffs: current block offset in the file.
   * @i: fid counting index.
   * @offset: current byte offset in the file.
   * @bytes: final read byte data.
   * @blk: currently read fid.
   */
  unsigned int ndata, noffs, i, offset;
  uint8_t *bytes, *blk;

  /* compute the number of data bytes to allocate. */
  ndata = nblk * szblk;
  *n = ndata;

  /* allocate an array for the bytes of data. */
  bytes = (uint8_t*) malloc(ndata * sizeof(uint8_t));

  /* check that the byte array was successfully allocated. */
  if (!bytes)
    return NULL;

  /* initialize the byte and block offsets. */
  offset = 0;
  noffs = 0;

  /* loop through the number of fids in the file. */
  for (i = 0; i < nblk; i++) {
    /* read the next fid from the file. */
    blk = bytes_read_block(fname, offset, szblk);

    /* check that the read was successful. */
    if (!blk) {
      /* nope. free the byte array and return failure. */
      free(bytes);
      return NULL;
    }

    /* move to the start of the next fid in the data file. all fids
     * in bruker ser files begin at 1.0 KiB boundaries.
     */
    offset += szblk;
    noffs = offset / 1024;
    while (offset % 1024)
      offset = ++noffs * 1024;

    /* copy the fid bytes into the byte array and free them. */
    memcpy(bytes + i * szblk, blk, szblk);
    free(blk);
  }

  /* return the read byte array. */
  return bytes;
}

/* bytes_read_varian(): reads fid byte data from a varian fid file.
 * @fname: the input filename.
 * @nblk: the number of fids in the file.
 * @szblk: the number of bytes per fid.
 * @offblk: the number of bytes per block header.
 * @offhead: the number of bytes in the file header.
 * @n: pointer to the output byte count.
 */
uint8_t *bytes_read_varian (const char *fname,
                            unsigned int nblk,
                            unsigned int szblk,
                            unsigned int offblk,
                            unsigned int offhead,
                            unsigned int *n) {
  /* declare a few required variables:
   * @ndata: number of bytes of real data in the file.
   * @i: fid counting index.
   * @offset: current byte offset in the file.
   * @bytes: final read byte data.
   * @blk: currently read fid.
   */
  int ndata, i, offset;
  uint8_t *bytes, *blk;

  /* compute the number of data bytes to allocate. */
  ndata = nblk * szblk;
  *n = ndata;

  /* allocate an array for the bytes of data. */
  bytes = (uint8_t*) malloc(ndata * sizeof(uint8_t));

  /* check that the byte array was successfully allocated. */
  if (!bytes)
    return 0;

  /* initialize the byte offset. */
  offset = offhead;

  /* loop through the number of blocks in the file. */
  for (i = 0; i < nblk; i++) {
    /* move past the block header bytes. */
    offset += offblk;

    /* read the next block from the file. */
    blk = bytes_read_block(fname, offset, szblk);

    /* check that the read was successful. */
    if (!blk) {
      /* nope. free the byte array and return failure. */
      free(bytes);
      return NULL;
    }

    /* move to the start of the next block in the data file.
     */
    offset += szblk;

    /* copy the fid bytes into the byte array and free them. */
    memcpy(bytes + i * szblk, blk, szblk);
    free(blk);
  }

  /* return the read byte array. */
  return bytes;
}

