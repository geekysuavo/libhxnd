
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

/* define a magic number to use when determining byte order of binary-format
 * NMR datum files. (in L.E. = 'HXNDDATA')
 */
#define NMR_DATUM_MAGIC  0x41544144444e5848

/* datum_type: enumerated type for datum raw byte data types.
 */
enum datum_type {
  DATUM_TYPE_UNDEFINED,
  DATUM_TYPE_BRUKER,
  DATUM_TYPE_VARIAN,
  DATUM_TYPE_PIPE,
  DATUM_TYPE_HXND,
  DATUM_TYPE_TEXT
};

/* datum_dim: single dimension of parameters for acquired NMR data.
 */
typedef struct {
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
   */
  unsigned int cx, nus, ft;

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

  /* @type: raw byte data type. */
  enum byteorder endian;
  enum datum_type type;

  /* @dims: array of per-dimension parameter values.
   * @nd: number of dimensions.
   */
  datum_dim *dims;
  unsigned int nd;

  /* @array: raw time-domain and/or frequency-domain NMR data.
   * @array_ok: whether the array has been allocated
   */
  int array_alloc;
  hx_array array;
}
datum;

/* include the nmr format headers. */
#include <hxnd/nmr-bruker.h>
#include <hxnd/nmr-varian.h>
#include <hxnd/nmr-pipe.h>

/* function declarations: */

void datum_init (datum *D);

void datum_free (datum *D);

enum datum_type datum_lookup_type (const char *name);

enum datum_type datum_guess_type (const char *fname);

int datum_print (datum *D, const char *fname);

int datum_check_magic (const char *fname);

int datum_fwrite_formatted (datum *D, FILE *fh, enum datum_type fmt);

int datum_fwrite_text (datum *D, FILE *fh);

int datum_fwrite (datum *D, FILE *fh);

int datum_fread (datum *D, FILE *fh, int read_array);

int datum_save (datum *D, const char *fname);

int datum_load (datum *D, const char *fname, int load_array);

int datum_get_dim_parameter (datum *D, const char *name, unsigned int d,
                             void *parm);

int datum_set_dim_parameter (datum *D, const char *name, unsigned int d,
                             const char *parm);

int datum_reorder_dims (datum *D, int *order);

int datum_refactor_array (datum *D);

int datum_read_array (datum *D);

int datum_free_array (datum *D);

int datum_resize_array (datum *D, int *sz);

int datum_slice_array (datum *D, int *lower, int *upper);

int datum_fill (datum *D, const char *fname);

#endif /* __HXND_NMR_DATUM_H__ */

