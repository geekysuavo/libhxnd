
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
#ifndef __HXND_NMR_VARIAN_H__
#define __HXND_NMR_VARIAN_H__

/* include the n-dimensional math header. */
#include <hxnd/hx.h>

/* include the byte-level data, nmr datum, and string headers. */
#include <hxnd/nmr-datum.h>
#include <hxnd/bytes.h>
#include <hxnd/str.h>

/* define constant parameter type characters for procpar file parsing.
 */
#define VARIAN_PARMTYPE_INT     'i'
#define VARIAN_PARMTYPE_INTS    'I'
#define VARIAN_PARMTYPE_FLOAT   'f'
#define VARIAN_PARMTYPE_STRING  's'

/* define bit masks for file header status words.
 */
#define VARIAN_HDR_S_DATA     0x0001
#define VARIAN_HDR_S_SPEC     0x0002
#define VARIAN_HDR_S_32       0x0004
#define VARIAN_HDR_S_FLOAT    0x0008
#define VARIAN_HDR_S_COMPLEX  0x0010
#define VARIAN_HDR_S_HYPERCX  0x0020
#define VARIAN_HDR_S_ACQPAR   0x0080
#define VARIAN_HDR_S_SECND    0x0100
#define VARIAN_HDR_S_TRANSF   0x0200
#define VARIAN_HDR_S_NP       0x0800
#define VARIAN_HDR_S_NF       0x1000
#define VARIAN_HDR_S_NI       0x2000
#define VARIAN_HDR_S_NI2      0x4000

/* varian_hdr_file: structure definition for varian file headers.
 */
struct varian_hdr_file {
  uint32_t nblocks;  /* number of data blocks. */
  uint32_t ntraces;  /* number of traces per block. */
  uint32_t np;       /* number of elements per trace. */
  uint32_t ebytes;   /* number of bytes per element. */
  uint32_t tbytes;   /* number of bytes per trace. */
  uint32_t bbytes;   /* number of bytes per block. */
  uint16_t vers_id;  /* software version, file_id status bits. */
  uint16_t status;   /* status of whole file. */
  uint32_t nheaders; /* number of block headers per block. */
};

/* varian_hdr_blk: structure definition for varian block headers.
 */
struct varian_hdr_blk {
  uint16_t scale;    /* scaling factor. */
  uint16_t status;   /* status of data in block. */
  uint16_t index;    /* block index. */
  uint16_t mode;     /* mode of data in block. */
  uint32_t ctcount;  /* ct value for fid. */
  float lpval;       /* f2 (2D-f1) left phase in phasefile. */
  float rpval;       /* f2 (2D-f1) right phase in phasefile. */
  float lvl;         /* level drift compensation. */
  float tlt;         /* tilt drift compensation. */
};

/* varian_hdr_ext: structure definition for varian extended block headers.
 */
struct varian_hdr_ext {
  uint16_t s_spare1;
  uint16_t status;    /* status word for block header. */
  uint16_t s_spare2;
  uint16_t s_spare3;
  uint32_t l_spare1;
  float lpval1;       /* 2D-f2 left phase. */
  float rpval1;       /* 2D-f2 right phase. */
  float f_spare1;
  float f_spare2;
};

/* function declarations: */

int varian_check_dir (const char *dname);

int varian_read_parms (const char *fname, unsigned int n, ...);

int varian_read_hdr_file (const char *fname, enum byteorder *endianness,
                          struct varian_hdr_file *hdr);

int varian_read (const char *fname, hx_array *x);

unsigned int varian_count_dims (const char *fname);

int varian_fill_datum (const char *dname, datum *D);

#endif /* __HXND_NMR_VARIAN_H__ */

