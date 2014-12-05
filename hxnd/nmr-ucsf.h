
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
#ifndef __HXND_NMR_UCSF_H__
#define __HXND_NMR_UCSF_H__

/* include the n-dimensional math header. */
#include <hxnd/hx.h>

/* include the byte-level data and nmr datum headers. */
#include <hxnd/nmr-datum.h>
#include <hxnd/bytes.h>

/* define a magic string used to check existence of ucsf files. */
#define UCSF_NUM_MAGIC  10
#define UCSF_MAGIC      "UCSF NMR"

/* ucsf_file_header: 180-byte file header of data contained in a ucsf-format
 * file.
 */
struct ucsf_file_header {
  /* (0..9) @ftype: ucsf header string.
   * (10) @ndims: dimensionality of spectrum.
   * (11) @ncomp: number of data components.
   * (12) @pad0
   * (13) @fmtver: format version number.
   */
  char ftype[UCSF_NUM_MAGIC];
  uint8_t ndims;
  uint8_t ncomp;
  uint8_t pad0;
  uint8_t fmtver;

  /* (14..179) @pad_end */
  uint8_t pad_end[166];
};

/* ucsf_dim_header: 128-byte dimension header of data contains in a ucsf
 * format file.
 */
struct ucsf_dim_header {
  /* (0..7) @nuc: nucleus name string.
   * (8..11) @npts: number of points.
   * (12..15) @pad0
   * (16..19) @sztile: tile size.
   * (20..23) @carrier: spectrometer frequency.
   * (24..27) @width: spectral width.
   * (28..31) @center: spectral center.
   */
  char nuc[8];
  uint32_t npts;
  uint32_t pad0;
  uint32_t sztile;
  float carrier;
  float width;
  float center;

  /* (32..127) @pad_end */
  uint8_t pad_end[96];
};

/* function declarations: */

int ucsf_check_magic (const char *fname);

int ucsf_read_header (const char *fname,
                      enum byteorder *endianness,
                      struct ucsf_file_header *hdr,
                      struct ucsf_dim_header **dims);

int ucsf_read (const char *fname, hx_array *x);

int ucsf_fill_datum (const char *fname, datum *D);

int ucsf_fwrite_datum (datum *D, FILE *fh);

#endif /* __HXND_NMR_UCSF_H__ */

