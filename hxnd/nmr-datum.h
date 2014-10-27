
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

/* datum_dim: single dimension of parameters for acquired NMR data.
 */
typedef struct {
  /* size parameters:
   * @sz: array of (actual) dimension sizes.
   * @td: array of (time-domain) dimension sizes.
   * @tdunif: array of (uniform time-domain) dimension sizes.
   */
  unsigned int sz, td, tdunif;

  /* frequency-domain parameters:
   */
  real carrier, width, offset;

  /* status flags:
   * @cx: complex (1) or real (0).
   * @nus: nonuniform (1) or uniform (0).
   */
  unsigned int cx, nus;

  /* @nuc: nucleus string. */
  char nuc[8];
}
datum_dim;

/* datum: data type for acquired NMR data.
 *
 * datum structures hold n-dimensional NMR datasets, and are essentially
 * encapsulations of hypercomplex multidimensional arrays with additional
 * metadata that describes all relevant features of the data.
 *
 * FIXME: complete the 'datum' structure type.
 */
typedef struct {
  /* @dims: array of per-dimension parameter values.
   * @nd: number of dimensions.
   */
  datum_dim *dims;
  unsigned int nd;

  /* @array: raw time-domain and/or frequency-domain NMR data. */
  hx_array array;
}
datum;

/* function declarations: */

int datum_print (datum *D, const char *fname);

int datum_save (datum *D, const char *fname);

int datum_load (datum *D, const char *fname);

int datum_reorder_dims (datum *D, int *order);

#endif /* __HXND_NMR_DATUM_H__ */

