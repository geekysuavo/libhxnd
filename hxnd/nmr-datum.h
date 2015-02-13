
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
#ifndef __HXND_NMR_DATUM_H__
#define __HXND_NMR_DATUM_H__

/* include the n-dimensional math header. */
#include <hxnd/hx.h>

/* include the byte-level data header. */
#include <hxnd/bytes.h>

/* define a value for datum core array dimension indices that indicates no
 * corresponding array dimension exists.
 */
#define DATUM_DIM_INVALID  (-1)

/* datum_type: enumerated type for datum raw byte data types.
 */
enum datum_type {
  DATUM_TYPE_UNDEFINED,
  DATUM_TYPE_HXND,
  DATUM_TYPE_TEXT,
  DATUM_TYPE_BRUKER,
  DATUM_TYPE_VARIAN,
  DATUM_TYPE_PIPE,
  DATUM_TYPE_UCSF,
  DATUM_TYPE_NV,
  DATUM_TYPE_RNMRTK
};

/* datum_dim: single dimension of parameters for acquired NMR data.
 */
typedef struct {
  /* array-datum correspondence parameters:
   * @d: array algebraic dimension, or -1 for none.
   * @k: array topological dimension. must be valid.
   */
  int d, k;

  /* size parameters:
   * @sz: array of (actual) dimension sizes.
   * @td: array of (time-domain) dimension sizes.
   * @tdunif: array of (uniform time-domain) dimension sizes.
   */
  unsigned int sz, td, tdunif;

  /* status flags:
   * @cx: complex (1) or real (0).
   * @nus: nonuniform (1) or uniform (0).
   * @ft: frequency (1) or time (0) domain.
   * @alt: sign alternation required (1) or not (0).
   * @neg: imaginary negation required (1) or not (0).
   * @genh: gradient-enhanced (1) or not (0).
   */
  unsigned int cx, nus, ft, alt, neg, genh;

  /* spectral parameters:
   */
  real carrier, width, offset;

  /* @nuc: nucleus string. */
  char nuc[8];
}
datum_dim;

/* datum: data type for acquired NMR data.
 *
 * datum structures hold n-dimensional NMR datasets, and are essentially
 * encapsulations of hypercomplex multidimensional arrays with additional
 * metadata that describes all relevant features of the data.
 */
typedef struct {
  /* @fname: raw byte data filename string. */
  char *fname;

  /* @endian: raw data byte ordering.
   * @type: file type for array reads.
   */
  enum byteorder endian;
  enum datum_type type;

  /* @epoch: time since the epoch of acquisition or datum creation.
   */
  time_t epoch;

  /* @dims: array of per-dimension parameter values.
   * @nd: number of dimensions.
   */
  datum_dim *dims;
  unsigned int nd;

  /* @sched: array of sampled indices generated from a nonuniform schedule.
   * @d_sched: number of values per sampled index in @sched.
   * @n_sched: number of sampled indices in @sched.
   */
  int d_sched, n_sched;
  int *sched;

  /* @grpdelay: bruker group delay correction value.
   */
  real grpdelay;

  /* @array: raw time-domain and/or frequency-domain NMR data.
   * @array_ok: whether the array has been allocated
   */
  int array_alloc;
  hx_array array;
}
datum;

/* include the nmr format headers. */
#include <hxnd/nmr-hxnd.h>
#include <hxnd/nmr-text.h>
#include <hxnd/nmr-bruker.h>
#include <hxnd/nmr-varian.h>
#include <hxnd/nmr-pipe.h>
#include <hxnd/nmr-ucsf.h>
#include <hxnd/nmr-nv.h>
#include <hxnd/nmr-rnmrtk.h>

/* function declarations (nmr-datum.c): */

int datum_read_sched (datum *D, const char *fname);

/* function declarations (nmr-datum-mem.c): */

void datum_init (datum *D);

void datum_free (datum *D);

/* function declarations (nmr-datum-io.c): */

int datum_load (datum *D, const char *fname);

int datum_print (datum *D, const char *fname);

/* function declarations (nmr-datum-dims.c): */

int datum_dims_getparm (datum *D, const char *name,
                        unsigned int d,
                        void *parm);

int datum_dims_setparm (datum *D, const char *name,
                        unsigned int d,
                        const char *parm);

int datum_dims_realloc (datum *D, unsigned int nd);

int datum_dims_reorder (datum *D, int *order);

/* function declarations (nmr-datum-array.c): */

int datum_array_refactor (datum *D);

int datum_array_alloc (datum *D);

int datum_array_read (datum *D);

int datum_array_free (datum *D);

int datum_array_resize (datum *D, int *sz);

int datum_array_slice (datum *D, int *lower, int *upper);

int datum_array_project (datum *D, int dim, hx_array_projector_cb projector);

/* function declarations (nmr-datum-type.c): */

const char *datum_type_name (enum datum_type typ);

const char *datum_type_desc (enum datum_type typ);

enum datum_type datum_type_lookup (const char *name);

enum datum_type datum_type_guess (const char *fname);

int datum_type_decode (datum *D, const char *fname, enum datum_type typ);

int datum_type_encode (datum *D, const char *fname, enum datum_type typ);

int datum_type_array (datum *D, enum datum_type typ);

int datum_type_post (datum *D, enum datum_type typ);

#endif /* __HXND_NMR_DATUM_H__ */

