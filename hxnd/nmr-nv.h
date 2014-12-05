
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

/* ensure once-only inclusion. */
#ifndef __HXND_NMR_NV_H__
#define __HXND_NMR_NV_H__

/* include the n-dimensional math header. */
#include <hxnd/hx.h>

/* include the byte-level data and nmr datum headers. */
#include <hxnd/nmr-datum.h>
#include <hxnd/bytes.h>

/* define a magic string used to check existence of nmrview files. */
#define NV_MAGIC      0

/* nv_dim_header: dimension header of data contains in an nmrview
 * format file.
 */
struct nv_dim_header {
  /* (0) @sz: FIXME
   * (1) @szblk: FIXME
   * (2) @nblk: FIXME
   * (3) @offblk: FIXME
   * (4) @maskblk: FIXME
   * (5) @ptoff: FIXME
   */
  int32_t sz;
  int32_t szblk;
  int32_t nblk;
  int32_t offblk;
  int32_t maskblk;
  int32_t ptoff;

  /* (6) @sf: FIXME
   * (7) @sw: FIXME
   * (8) @refpt: FIXME
   * (9) @ref: FIXME
   * (10) @refunits: FIXME
   */
  float sf;
  float sw;
  float refpt;
  float ref;
  int32_t refunits;

  /* (11) @foldup: FIXME
   * (12) @folddown: FIXME
   */
  float foldup;
  float folddown;

  /* (13..28) @label: dimension label string. */
  char label[16];

  /* (29..43) @pad_end */
  int32_t pad_end[15];
};

/* nv_file_header: complete file header of data contained in an nmrview
 * format file.
 */
struct nv_file_header {
  /* (0) @magic: magic number.
   * (1..2) @pad0
   * (3) @fhdrsz: FIXME
   * (4) @bhdrsz: FIXME
   * (5) @blkelem: FIXME
   * (6) @ndims: FIXME
   * (7) @temp: collection temperature.
   * (8..39) @sequence: pulse sequence string.
   * (40..199) @comment: comment string.
   * (200) @month: collection month.
   * (201) @day: collection day of month.
   * (202) @year: collection year.
   */
  int32_t magic;
  int32_t pad0[2];
  int32_t fhdrsz;
  int32_t bhdrsz;
  int32_t blkelem;
  int32_t ndims;
  float temp;
  char sequence[32];
  char comment[160];
  int32_t month;
  int32_t day;
  int32_t year;

  /* (203..399) @pad1 */
  int32_t pad1[197];

  /* (400..end) @dims: dimension fields. */
  struct nv_dim_header dims[8];
};

/* function declarations: */

int nv_check_magic (const char *fname);

int nv_read_header (const char *fname,
                    enum byteorder *endianness,
                    struct nv_file_header *hdr);

int nv_read (const char *fname, hx_array *x);

int nv_fill_datum (const char *fname, datum *D);

int nv_fwrite_datum (datum *D, FILE *fh);

#endif /* __HXND_NMR_NV_H__ */

