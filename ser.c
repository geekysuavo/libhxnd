
/* ndmath: A framework for loading free induction decay series data.
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

/* include the byte loading and series headers. */
#include "bytes.h"
#include "ser.h"

/* ser_word(): converts a single word into the correct format for inclusion
 * into a hypercomplex array structure.
 * @bytes: the bytes of the word.
 * @sz: the byte count of the word.
 * @flt: if the word is floating point.
 */
real ser_word (uint8_t *bytes, int sz, int flt) {
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

/* ser_toarray(): converts bytes loaded from a bruker/varian fid/ser
 * serial file into a one-dimensional real array.
 * @bytes: the byte array to convert into an array.
 * @nbytes: the number of bytes in the array.
 * @endianness: the endianness of the data.
 * @wordsz: number of bytes per data word.
 * @flt: whether each word is a float (1) or int (0).
 * @x: the final output array.
 */
int ser_toarray (uint8_t *bytes, unsigned int nbytes,
                 int endianness, int wordsz, int flt,
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
  if (endianness == BYTES_ENDIAN_BIG && wordsz > 1) {
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
    x->x[k] = ser_word(bytes + i, wordsz, flt);
  }

  /* return success. */
  return 1;
}

/* ser_deinterlace(): convert the largest dimension (the last index) of an
 * array to the next level of algebraic dimensionality (i.e. x->d++). for
 * example, perform:
 *
 *   (d=0, k=1, sz={2*n})       -> (d=1, k=1, sz={n})
 *   (d=1, k=2, sz={n, 2*m})    -> (d=2, k=2, sz={n, m})
 *   (d=2, k=3, sz={n, m, 2*l}) -> (d=3, k=3, sz={n, m, l})
 *
 * interlacing (d = 0 -> 1):
 *   [1, u1, ...]
 *    -> [(1, u1), ...]
 *
 * interlacing (d = 1 -> 2):
 *   [(1, u1), (u2, u1u2), ...]
 *    -> [(1, u1, u2, u1u2), ...]
 *
 * interlacing (d = 2 -> 3):
 *   [(1, u1, u2, u1u2), (u3, u1u3, u2u3, u1u2u3), ...]
 *    -> [(1, u1, u2, u1u2, u3, u1u3, u2u3, u1u2u3), ...]
 *
 * it is assumed that each alternating point in the topmost array dimension
 * contains the next set of coefficients for the extra algebraic dimensions
 * created by promotion to the next greatest dimensionality.
 */
int ser_deinterlace (hx_array *x) {
  /* declare a few required variables:
   * @ktop: the (topmost) array dimension to act upon.
   * @dnew: the new algebraic dimensionality of the array.
   * @i: a general-purpose dim/coefficient loop counter.
   * @di: coefficient index offset during reshuffling.
   * @arri: the input index array.
   * @arro: the output index array.
   * @sznew: the new size array.
   * @idxi: the input linear index.
   * @idxo: the output linear index.
   */
  int ktop, dnew, i, di, *arri, *arro, *sznew, idxi, idxo;

  /* get the topmost array dimension. */
  ktop = x->k - 1;

  /* check that the topmost array dimension is valid. */
  if (ktop < 0)
    return 0;

  /* check that the topmost array dimension size is divisible by two. */
  if (x->sz[ktop] % 2)
    return 0;

  /* allocate a new size array. */
  sznew = hx_array_index_alloc(x->k);

  /* check that allocation was successful. */
  if (!sznew)
    return 0;

  /* copy the current array sizes into the new size array. */
  memcpy(sznew, x->sz, x->k * sizeof(int));

  /* get the new dimensionality of the array. */
  dnew = x->d + 1;

  /* attempt to promote the array from real to complex. */
  if (!hx_array_resize(x, dnew, x->k, sznew))
    return 0;

  /* initialize a set of multidimensional index arrays. */
  idxi = idxo = 0;
  arri = hx_array_index_alloc(x->k);
  arro = hx_array_index_alloc(x->k);
  if (!arri || !arro)
    return 0;

  /* iterate over the points of the array. */
  do {
    /* copy the relevant indices of the index array. */
    for (i = 0; i < ktop; i++)
      arro[i] = arri[i];

    /* determine which set of coefficients we're on (i.e. 'odd' or 'even').
     */
    if (arri[ktop] % 2 == 0) {
      /* even. this will copy the first half of the coefficients. */
      arro[ktop] = arri[ktop] / 2;
      di = 0;
    }
    else {
      /* odd. this will copy the second half of the coefficients. */
      arro[ktop] = (arri[ktop] - 1) / 2;
      di = x->n / 2;
    }

    /* pack the destination index array into a linear index. */
    hx_array_index_pack(x->k, sznew, arro, &idxo);

    /* loop over the first half of the coefficients, copying their values
     * into their final destinations.
     */
    for (i = 0; i < x->n / 2; i++)
      x->x[i + di + x->n * idxo] = x->x[i + x->n * idxi];

    /* increment the input linear index. */
    idxi++;
  } while (hx_array_index_inc(x->k, sznew, &arri));

  /* free the array indices. */
  free(arri);
  free(arro);

  /* shrink the topmost array dimension by half. */
  sznew[ktop] /= 2;

  /* shrink the geometry of the array. */
  if (!hx_array_resize(x, dnew, x->k, sznew))
    return 0;

  /* free the new size array. */
  free(sznew);

  /* return success. */
  return 1;
}

