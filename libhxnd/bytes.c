
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

/* bytes_conv_real_u64_t: union definition for converting between a real
 * floating point value and its binary representation.
 */
typedef union {
  real fval;
  uint64_t ival;
}
bytes_conv_real_u64_t;

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

/* bytes_get_native(): returns the native byte ordering.
 */
enum byteorder bytes_get_native (void) {
  /* ensure that the native endianness has been determined. */
  bytes_init();

  /* return the native endianness. */
  return bytes_native_endianness;
}

/* bytes_get_nonnative(): returns the opposite of native byte ordering.
 */
enum byteorder bytes_get_nonnative (void) {
  /* ensure that the native endianness has been determined. */
  bytes_init();

  /* return the non-native endianness. */
  if (bytes_native_endianness == BYTES_ENDIAN_BIG)
    return BYTES_ENDIAN_LITTLE;
  else
    return BYTES_ENDIAN_BIG;
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

/* bytes_swap(): swaps the bytes of an arbitrarily-sized array of words.
 * @bytes: array of bytes to reorder.
 * @n: number of words in the array.
 * @sz: size of each word, in bytes.
 */
void bytes_swap (uint8_t *bytes, unsigned int n, unsigned int sz) {
  /* declare a few required variables:
   * @i, @j: loop counters.
   * @jmax: the @j-loop boundary.
   * @swp: the swapped byte value.
   */
  unsigned int i, j, jmax, nbytes;
  uint8_t swp;

  /* compute the inner loop boundary. this is just a really slick way of
   * determining how many byte pairs need to be exchanged, given a word
   * size.
   */
  jmax = ((sz + 1) & ~0x01) / 2;

  /* compute the total number of bytes in the array. */
  nbytes = n * sz;

  /* loop over the words, in byte units. */
  for (i = 0; i < nbytes; i += sz) {
    /* loop over the bytes of the word. */
    for (j = 0; j < jmax; j++) {
      /* swap bytes. */
      swp = bytes[i + j];
      bytes[i + j] = bytes[i + sz - j - 1];
      bytes[i + sz - j - 1] = swp;
    }
  }
}

/* bytes_fexist(): return whether the regular file @fname exists.
 * @fname: the filename to check.
 */
int bytes_fexist (const char *fname) {
  /* declare a required variable. */
  int ret = 0;
  FILE *fh;

  /* open the input file. */
  fh = fopen(fname, "rb");

  /* check if the file was opened. */
  if (fh) {
    /* close the file and raise the return value. */
    fclose(fh);
    ret = 1;
  }
  else {
    /* clear the error number. */
    traceback_clear();
  }

  /* return the result. */
  return ret;
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
    throw("failed to open '%s'", fname);

  /* move to the end of the file. */
  if (fseek(fh, 0, SEEK_END))
    throw("failed to seek '%s'", fname);

  /* read the file size. */
  n = (unsigned int) ftell(fh);

  /* move back to the beginning of the file. */
  if (fseek(fh, 0, SEEK_SET))
    throw("failed to seek '%s'", fname);

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
  if (!fh) {
    /* raise an error and return nothing. */
    raise("failed to open '%s'", fname);
    return NULL;
  }

  /* allocate memory for the read bytes. */
  bytes = (uint8_t*) malloc(n * sizeof(uint8_t));

  /* check that the memory was allocated. */
  if (!bytes) {
    /* raise an error and return nothing. */
    raise("failed to allocate %u bytes", n);
    return NULL;
  }

  /* move back to the beginning of the file. */
  if (fseek(fh, offset, SEEK_SET)) {
    /* free the byte array and return nothing. */
    raise("failed to seek '%s'", fname);
    free(bytes);
    return NULL;
  }

  /* read the bytes in from the input file. */
  if (fread(bytes, sizeof(uint8_t), n, fh) != n) {
    /* free the byte array and return nothing. */
    raise("failed to read %u bytes from '%s'", n, fname);
    free(bytes);
    return NULL;
  }

  /* close the input file and return the read data. */
  fclose(fh);
  return bytes;
}

/* bytes_real_to_u64(): convert a real floating point value to a u64.
 */
uint64_t bytes_real_to_u64 (const real x) {
  /* declare a conversion union. */
  bytes_conv_real_u64_t conv;

  /* store the input value into the union. */
  conv.fval = x;

  /* return the output value. */
  return conv.ival;
}

/* bytes_u64_to_real(): convert a u64 to a real floating point value.
 */
real bytes_u64_to_real (const uint64_t x) {
  /* declare a conversion union. */
  bytes_conv_real_u64_t conv;

  /* store the input value into the union. */
  conv.ival = x;

  /* return the output value. */
  return conv.fval;
}

/* bytes_unpack(): convert a single raw data word into the native hx format.
 * @bytes: the bytes of the word.
 * @sz: the byte count of the word.
 * @isflt: if the word is floating point.
 */
real bytes_unpack (uint8_t *bytes, int sz, int isflt) {
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
  if (isflt) {
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
      x = (real) i8;
    }
    else if (sz == sizeof(int16_t)) {
      /* build a two-byte word. */
      i16 = (bytes[1] << 8) |
             bytes[0];

      /* convert to float. */
      x = (real) i16;
    }
    else if (sz == sizeof(int32_t)) {
      /* build a four-byte word. */
      i32 = (bytes[3] << 24) |
            (bytes[2] << 16) |
            (bytes[1] << 8) |
             bytes[0];

      /* convert to float. */
      x = (real) i32;
    }
    else if (sz == sizeof(int64_t)) {
      /* build an eight-byte float word. */
      i64 = ((int64_t) bytes[7] << 56) | ((int64_t) bytes[6] << 48) |
            ((int64_t) bytes[5] << 40) | ((int64_t) bytes[4] << 32) |
            (bytes[3] << 24) | (bytes[2] << 16) |
            (bytes[1] << 8) |
             bytes[0];

      /* convert to the final value. */
      x = (real) i64;
    }
  }

  /* return the converted real value. */
  return x;
}

/* bytes_pack(): build a single raw data word from the native hx format.
 * @value: the input native hx-format value.
 * @bytes: the generated raw bytes of the word.
 * @sz: the byte count of the output word.
 * @isflt: if the output word is floating point.
 */
int bytes_pack (real value, uint8_t *bytes, int sz, int isflt) {
  /* declare variables for packing.
   */
  int16_t i16;
  int32_t i32;
  int64_t i64;
  float f32;
  double f64;

  /* convert based on whether we are in integer or float format. */
  if (isflt) {
    /* check the size of the float. */
    if (sz == sizeof(float)) {
      /* convert to a four-byte float word. */
      f32 = (float) value;
      memcpy(bytes, &f32, sz);
    }
    else if (sz == sizeof(double)) {
      /* convert to an eight-byte float word. */
      f64 = (double) value;
      memcpy(bytes, &f64, sz);
    }
    else
      throw("cannot pack %d-byte floats", sz);
  }
  else {
    /* check the size of the integer. */
    if (sz == sizeof(int8_t)) {
      /* convert to a one-byte integer. */
      bytes[0] = (int8_t) value;
    }
    else if (sz == sizeof(int16_t)) {
      /* convert to a two-byte integer. */
      i16 = (int16_t) value;
      bytes[0] = BYTES_BYTE0(i16);
      bytes[1] = BYTES_BYTE1(i16);
    }
    else if (sz == sizeof(int32_t)) {
      /* convert to a four-byte integer. */
      i32 = (int32_t) value;
      bytes[0] = BYTES_BYTE0(i32);
      bytes[1] = BYTES_BYTE1(i32);
      bytes[2] = BYTES_BYTE2(i32);
      bytes[3] = BYTES_BYTE3(i32);
    }
    else if (sz == sizeof(int64_t)) {
      /* convert to an eight-byte integer. */
      i64 = (int64_t) value;
      i32 = (int32_t) BYTES_RSHIFT4(i64);

      /* build the first four bytes. */
      bytes[0] = BYTES_BYTE0(i64);
      bytes[1] = BYTES_BYTE1(i64);
      bytes[2] = BYTES_BYTE2(i64);
      bytes[3] = BYTES_BYTE3(i64);

      /* build the next four bytes. */
      i32 = (int32_t) BYTES_RSHIFT4(i64);
      bytes[4] = BYTES_BYTE0(i32);
      bytes[5] = BYTES_BYTE1(i32);
      bytes[6] = BYTES_BYTE2(i32);
      bytes[7] = BYTES_BYTE3(i32);
    }
    else
      throw("cannot pack %d-byte integers", sz);
  }

  /* return success. */
  return 1;
}

