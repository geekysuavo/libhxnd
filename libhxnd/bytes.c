
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

/* include the byte-level data header. */
#include <hxnd/bytes.h>

/* bytes_native_endianness: a private variable in the bytes_* namespace that
 * identifies the native byte ordering of words on the current machine.
 */
enum byteorder bytes_native_endianness = BYTES_ENDIAN_AUTO;

/* bytes_init(): initializes the native endianness variable that may be used
 * later to determine if byte swaps are required.
 */
void bytes_init (void) {
  /* check if the endianness has already been determined. */
  if (bytes_native_endianness != BYTES_ENDIAN_AUTO)
    return;

  /* declare a unioned two-byte word and set it's value. */
  union {
    uint16_t ival;
    uint8_t bytes[2];
  } u = {0x0100};

  /* check if the first byte of the word is nonzero. */
  if (u.bytes[0] == 1) {
    /* the machine is using big-endian. */
    bytes_native_endianness = BYTES_ENDIAN_BIG;
  }
  else {
    /* the machine is using little-endian. */
    bytes_native_endianness = BYTES_ENDIAN_LITTLE;
  }
}

/* bytes_native(): returns whether the specified byte ordering is native
 * to the current machine.
 */
int bytes_native (enum byteorder endianness) {
  /* ensure that the native endianness has been determined. */
  bytes_init();

  /* check the value against the native value. */
  return (endianness == bytes_native_endianness);
}

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

/* bytes_toword(): converts a single word into the correct format for
 * inclusion into a hypercomplex array structure.
 * @bytes: the bytes of the word.
 * @sz: the byte count of the word.
 * @flt: if the word is floating point.
 */
real bytes_toword (uint8_t *bytes, int sz, int flt) {
  /* declare required variables to ensure proper signed integer support.
   * if the specified word size is either that of 'int' or 'short', then
   * each data word is down-cast to that type and then re-cast to a 'long',
   * which is a quick and dirty way to support signed integers of all
   * three sizes.
   */
  int8_t i8;
  int16_t i16;
  int32_t i32;
  int64_t i64;

  /* declare required variables to ensure proper float support.
   */
  real x;
  float f32;
  double f64;

  /* declare required variables to ensure proper half-precision support.
   */
  real sign;
  uint8_t expo;
  uint16_t scnd;

  /* initialize the converted value. */
  x = 0.0;

  /* convert based on whether we are in integer or float format. */
  if (flt) {
    /* check the size of the float. */
    if (sz == 2) {
      /* build a 16-bit word. */
      i16 = (bytes[1] << 8) | bytes[0];

      /* compute the sign, exponent and significand. */
      sign = (i16 & 0x8000 ? -1.0 : 1.0);
      expo = (i16 & 0x7c00) >> 10;
      scnd = (i16 & 0x03ff);

      /* determine how to interpret the value. this does not support
       * infinity or not-a-number.
       */
      if (expo == 0 && scnd) {
        /* compute a subnormal number. */
        x = sign * pow(2.0, -14.0);
        x *= ((real) scnd / 1024.0);
      }
      else {
        /* compute a normalized number. */
        x = sign * pow(2.0, (real) expo - 15.0);
        x *= (1.0 + (real) scnd / 1024.0);
      }
    }
    else if (sz == sizeof(float)) {
      /* build a four-byte float word and convert to real. */
      f32 = *(float*) bytes;
      x = (real) f32;
    }
    else if (sz == sizeof(double)) {
      /* build an eight-byte float word and convert to real. */
      f64 = *(double*) bytes;
      x = (real) f64;
    }
  }
  else {
    /* check the size of the integer. */
    if (sz == sizeof(int8_t)) {
      /* build a single-byte word and convert to float. */
      i8 = bytes[0];
      x = (real) i8 / 128.0;
    }
    if (sz == sizeof(int16_t)) {
      /* build a two-byte word. */
      i16 = (bytes[1] << 8) |
             bytes[0];

      /* convert to float. */
      x = (real) i16 / 32768.0;
    }
    else if (sz == sizeof(int32_t)) {
      /* build a four-byte word. */
      i32 = (bytes[3] << 24) |
            (bytes[2] << 16) |
            (bytes[1] << 8) |
             bytes[0];

      /* convert to float. */
      x = (real) i32 / 2147483648.0;
    }
    else if (sz == sizeof(int64_t)) {
      /* build an eight-byte float word. */
      i64 = ((int64_t) bytes[7] << 56) | ((int64_t) bytes[6] << 48) |
            ((int64_t) bytes[5] << 40) | ((int64_t) bytes[4] << 32) |
            (bytes[3] << 24) | (bytes[2] << 16) |
            (bytes[1] << 8) |
             bytes[0];

      /* convert to the final value. */
      x = (real) i64 / 9223372036854775808.0;
    }
  }

  /* return the converted real value. */
  return x;
}

/* bytes_toarray(): converts bytes loaded from a bruker/varian fid/ser
 * serial file into a one-dimensional real array.
 * @bytes: the byte array to convert into an array.
 * @nbytes: the number of bytes in the array.
 * @endianness: the endianness of the data.
 * @wordsz: number of bytes per data word.
 * @flt: whether each word is a float (1) or int (0).
 * @x: the final output array.
 */
int bytes_toarray (uint8_t *bytes, unsigned int nbytes,
                   enum byteorder endianness,
                   int wordsz, int flt,
                   hx_array *x) {
  /* declare a few required variables. */
  int nwords, topo[1], i, j, k;
  uint8_t swp;

  /* read the number of bytes in the input file and compute the number
   * of words in the serial file and words in the array.
   */
  nwords = nbytes / wordsz;
  topo[0] = nwords;

  /* check if byte swaps are required. */
  if (!bytes_native(endianness) && wordsz > 1) {
    /* loop over the byte array, per-word. */
    for (i = 0; i < nbytes; i += wordsz) {
      /* loop over the bytes of the word. */
      for (j = 0; j < wordsz / 2; j++) {
        /* swap bytes. */
        swp = bytes[i + j];
        bytes[i + j] = bytes[i + wordsz - j - 1];
        bytes[i + wordsz - j - 1] = swp;
      }
    }
  }

  /* allocate memory for a linear, real output array. */
  if (!hx_array_alloc(x, 0, 1, topo))
    return 0;

  /* copy the read, properly ordered byte data into the output array. */
  for (i = 0, k = 0; i < nbytes; i += wordsz, k++) {
    /* convert the data word into floating point. */
    x->x[k] = bytes_toword(bytes + i, wordsz, flt);
  }

  /* return success. */
  return 1;
}

