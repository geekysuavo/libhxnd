
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
#ifndef __HXND_NMR_RNMRTK_H__
#define __HXND_NMR_RNMRTK_H__

/* include the n-dimensional math header. */
#include <hxnd/hx.h>

/* include the byte-level data, string, and nmr datum headers. */
#include <hxnd/nmr-datum.h>
#include <hxnd/bytes.h>
#include <hxnd/str.h>

/* define the maximum number of dimensions supported by rnmrtk files. */
#define RNMRTK_MAXDIM   4
#define RNMRTK_MAXSUB  10

/* rnmrtk_quad: enumerated type describing the quadrature mode of a given
 * dimension in an RNMRTK data file.
 */
enum rnmrtk_quad {
  RNMRTK_QUAD_UNKNOWN    = 0x00,
  RNMRTK_QUAD_TPPI       = 0x01,
  RNMRTK_QUAD_STATES     = 0x02,
  RNMRTK_QUAD_STATESTPPI = 0x03
};

/* rnmrtk_parms: structure containing parsed parameters that correspond to
 * an RNMRTK data file.
 */
struct rnmrtk_parms {
  /* arguments parsed from parfile 'format' lines:
   * @endian: raw data byte ordering, little or big.
   * @isflt: whether the data is floating-point or integer.
   * @nheader: number of bytes in the data file header.
   * @reclen: number of real points per record.
   * @nbegin: number of bytes padding each data record.
   * @nend: number of bytes padding the end of the file.
   */
  enum byteorder endian;
  unsigned int isflt;
  unsigned int nheader;
  unsigned int reclen;
  unsigned int nbegin;
  unsigned int nend;

  /* arguments parsed from parfile 'dom' lines:
   * @ord: dimension ordering array, one-based.
   */
  int ord[RNMRTK_MAXDIM];

  /* arguments parsed from parfile 'n' lines:
   * @sz: array of sizes for each dimension.
   * @cx: array of complex flags for each dimension.
   */
  int sz[RNMRTK_MAXDIM], cx[RNMRTK_MAXDIM];

  /* arguments parsed from parfile 'layout' lines:
   * @layout: matrix of dimension,subdimension sizes.
   */
  int layout[RNMRTK_MAXDIM][RNMRTK_MAXSUB];

  /* arguments parsed from optional parfile lines:
   * @sf: spectrometer carrier frequencies.
   * @ppm: carrier offset values, in ppm.
   * @sw: spectral widths.
   */
  float sf[RNMRTK_MAXDIM];
  float ppm[RNMRTK_MAXDIM];
  float sw[RNMRTK_MAXDIM];
  enum rnmrtk_quad quad[RNMRTK_MAXDIM];
};

/* function declarations: */

int rnmrtk_check_file (const char *fname);

int rnmrtk_read_parms (const char *fname, struct rnmrtk_parms *par);

int rnmrtk_read (const char *fname, hx_array *x);

int rnmrtk_fill_datum (const char *fname, datum *D);

#endif /* __HXND_NMR_RNMRTK_H__ */

