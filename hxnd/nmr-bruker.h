
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
#ifndef __HXND_NMR_BRUKER_H__
#define __HXND_NMR_BRUKER_H__

/* include the n-dimensional math header. */
#include <hxnd/hx.h>

/* include the byte-level data, string, and nmr datum headers. */
#include <hxnd/nmr-datum.h>
#include <hxnd/bytes.h>
#include <hxnd/str.h>

/* define constant parameter type characters for acqus file parsing.
 */
#define BRUKER_PARMTYPE_INT     'i'
#define BRUKER_PARMTYPE_FLOAT   'f'
#define BRUKER_PARMTYPE_STRING  's'

/* bruker_aqmod: enumerated type for direct acquisition modes.
 */
enum bruker_aqmod {
  BRUKER_AQMOD_QF,    /* single-channel (real) detection. */
  BRUKER_AQMOD_QSIM,  /* simultaneous quadrature detection. */
  BRUKER_AQMOD_QSEQ,  /* sequential quadrature detection. */
  BRUKER_AQMOD_DQD    /* digital quadrature detection. */
};

/* bruker_fnmode: enumerated type for indirect acquisition modes.
 */
enum bruker_fnmode {
  BRUKER_FNMODE_UNDEFINED,
  BRUKER_FNMODE_QF,
  BRUKER_FNMODE_QSEQ,
  BRUKER_FNMODE_TPPI,
  BRUKER_FNMODE_STATES,
  BRUKER_FNMODE_STATESTPPI,
  BRUKER_FNMODE_GRADIENT
};

/* bruker_aqseq: (3+)-dimensional acquisition ordering scheme.
 */
enum bruker_aqseq {
  BRUKER_AQSEQ_321,
  BRUKER_AQSEQ_312
};

/* function declarations: */

int bruker_read_parms (const char *fname, unsigned int n, ...);

int bruker_read (const char *fname, enum byteorder endianness,
                 unsigned int nblk, unsigned int szblk,
                 hx_array *x);

int bruker_datum (const char *dname, datum *D);

#endif /* __HXND_NMR_BRUKER_H__ */

