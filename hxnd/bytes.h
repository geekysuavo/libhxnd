
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

/* ensure once-only inclusion. */
#ifndef __HXND_BYTES_H__
#define __HXND_BYTES_H__

/* include the n-dimensional math header. */
#include <hxnd/hx.h>

/* define macros to extract bytes from words and double-words.
 */
#define BYTES_BYTE0(x)     ((x) & 0x000000ff)
#define BYTES_BYTE1(x)    (((x) & 0x0000ff00) >> 8)
#define BYTES_BYTE2(x)    (((x) & 0x00ff0000) >> 16)
#define BYTES_BYTE3(x)    (((x) & 0xff000000) >> 24)

/* define macros to shift bytes by four.
 */
#define BYTES_LSHIFT4(x)  ((uint64_t) (x) << 32)
#define BYTES_RSHIFT4(x)  ((uint64_t) (x) >> 32)

/* define a macro to build a 16-bit word from bytes.
 */
#define BYTES_U16(b0, b1) \
  (((b1) << 8) | (b0))

/* define a macro to build a 32-bit word from bytes.
 */
#define BYTES_U32(b0, b1, b2, b3) \
  (((b3) << 24) | ((b2) << 16) | ((b1) << 8) | (b0))

/* define a macro to build a 64-bit word from bytes.
 */
#define BYTES_U64(b0, b1, b2, b3, b4, b5, b6, b7) \
  (BYTES_LSHIFT4(((b7) << 24) | ((b6) << 16) | ((b5) << 8) | (b4))| \
                (((b3) << 24) | ((b2) << 16) | ((b1) << 8) | (b0)))

/* define constants to specify the endianness of loaded serial files.
 */
#define BYTES_ENDIAN_AUTO    0x00
#define BYTES_ENDIAN_LITTLE  0x01
#define BYTES_ENDIAN_BIG     0x02

/* function declarations: */

void bytes_swap_u16 (uint16_t *x);

void bytes_swap_u32 (uint32_t *x);

void bytes_swap_u64 (uint64_t *x);

unsigned int bytes_size (const char *fname);

uint8_t *bytes_read_block (const char *fname,
                           unsigned int offset,
                           unsigned int n);

uint8_t *bytes_read_bruker (const char *fname,
                            unsigned int nblk,
                            unsigned int szblk,
                            unsigned int *n);

uint8_t *bytes_read_varian (const char *fname,
                            unsigned int nblk,
                            unsigned int szblk,
                            unsigned int offblk,
                            unsigned int offhead,
                            unsigned int *n);

int bytes_toarray (uint8_t *bytes, unsigned int nbytes,
                   int endianness, int wordsz, int flt,
                   hx_array *x);

#endif /* __HXND_BYTES_H__ */

