
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
#ifndef __HXND_HX_ARRAY_H__
#define __HXND_HX_ARRAY_H__

/* include the byte-level data header. */
#include <hxnd/bytes.h>

/* define a magic number to use when determining byte order of binary-format
 * hypercomplex multidimensional array files. (in L.E. = 'HXNDARRY')
 */
#define HX_ARRAY_MAGIC  0x59525241444e5848

/* define constants for forward (slice) and reverse (store) slicer operations.
 */
#define HX_ARRAY_SLICER_SLICE 0
#define HX_ARRAY_SLICER_STORE 1

/* define constants for forward (linearize) and reverse (tileize) tiler
 * operations.
 */
#define HX_ARRAY_TILER_FORWARD  0
#define HX_ARRAY_TILER_REVERSE  1

/* define constants for normal or reverse incrementation during tiling.
 */
#define HX_ARRAY_INCR_NORMAL   0
#define HX_ARRAY_INCR_REVERSE  1

/* hx_array: data type for nD arrays of hypercomplex nD numbers.
 *
 * all the above information still applies to arrays, but an extra layer of
 * complication is added to arrays, namely: they are arrays. instead of just
 * using literal arrays of hx_scalar structures, the array data is squashed
 * together into a single C array of real scalars.
 *
 * as it turns out, this squashing not only aids in memory compactness of
 * multidimensional datasets, but also in communicating array data between
 * the opencl host and device layers.
 */
typedef struct {
  /* d: dimensionality of the hypercomplex space.
   * n: number of coefficients per hypercomplex value (2**d).
   */
  int d, n;

  /* k: dimensionality of the array.
   * len: number of total array coefficients.
   * sz: array sizes along each array dimension.
   */
  int k, len, *sz;

  /* real coefficients. */
  real *x;

  /* multiplication table. this is also a shared table. */
  hx_algebra tbl;
}
hx_array;

/* hx_array_foreach_cb: callback function prototype for per-vector and
 * per-matrix operations that act on hypercomplex arrays.
 * @x: the host array pointer.
 * @y: the current slice from the host array.
 * @arr: the set of array indices for the first vector point.
 * @idx: the linear array index for the first vector point.
 * @vl: custom arguments, stored in a variable arguments list.
 */
typedef int (*hx_array_foreach_cb) (hx_array *x, hx_array *y,
                                    int *arr, int idx,
                                    va_list *vl);

/* hx_array_projector_cb: callback function prototype for array 'projection'
 * operations that act on hypercomplex arrays.
 * @y: the current slice from the host array.
 * @val: the result of the 'projection'.
 */
typedef int (*hx_array_projector_cb) (hx_array *y, real *val);

/* function declarations (hx-array.c): */

int hx_array_set_coeff (hx_array *x, int di, real value, ...);

int hx_array_complexify (hx_array *x, int genh);

int hx_array_real (hx_array *x, int d);

int hx_array_reshape (hx_array *x, int k, int *sz);

int hx_array_repack (hx_array *x, int ndiv);

int hx_array_shift (hx_array *x, int k, int amount);

/* function declarations (hx-array-mem.c): */

void hx_array_init (hx_array *x);

int hx_array_alloc (hx_array *x, int d, int k, int *sz);

int hx_array_copy (hx_array *dst, hx_array *src);

int hx_array_copy_real (hx_array *dst, hx_array *src);

void hx_array_free (hx_array *x);

/* function declarations (hx-array-io.c): */

int hx_array_print (hx_array *x, const char *fname);

int hx_array_check_magic (const char *fname);

int hx_array_fwrite (hx_array *x, FILE *fh);

int hx_array_fread (hx_array *x, FILE *fh);

int hx_array_save (hx_array *x, const char *fname);

int hx_array_load (hx_array *x, const char *fname);

/* function declarations (hx-array-rawio.c): */

int hx_array_fread_raw (FILE *fh, hx_array *x, enum byteorder endian,
                        unsigned int wordsz, unsigned int isflt,
                        unsigned int offhead, unsigned int offblk,
                        unsigned int nblks, unsigned int nwords,
                        unsigned int nalign);

int hx_array_fwrite_raw (FILE *fh, hx_array *x, enum byteorder endian,
                         unsigned int wordsz, unsigned int isflt);

/* function declarations (hx-array-topo.c): */

int hx_array_nnzdims (hx_array *x);

int hx_array_is_vector (hx_array *x);

int hx_array_is_matrix (hx_array *x);

int hx_array_is_cube (hx_array *x);

#define hx_array_assert_vector(x) \
  if (!hx_array_is_vector(x)) \
    throw("expected vector configuration for array %s", #x);

#define hx_array_assert_matrix(x) \
  if (!hx_array_is_matrix(x)) \
    throw("expected matrix configuration for array %s", #x);

#define hx_vector_len(x) \
  (hx_array_is_vector(x) ? x->sz[0] : 0)

#define hx_matrix_rows(x) \
  (hx_array_is_matrix(x) ? x->sz[0] : 0)

#define hx_matrix_cols(x) \
  (hx_array_is_matrix(x) ? x->sz[1] : 0)

int hx_array_compact (hx_array *x);

/* function declarations (hx-array-resize.c): */

int hx_array_is_real (hx_array *x);

int hx_array_resize (hx_array *x, int d, int k, int *sz);

/* function declarations (hx-array-slice.c): */

int hx_array_slicer (hx_array *x, hx_array *y,
                     int *lower, int *upper,
                     int dir);

#define hx_array_slice(x, y, l, u) \
  hx_array_slicer(x, y, l, u, HX_ARRAY_SLICER_SLICE)

#define hx_array_store(x, y, l, u) \
  hx_array_slicer(x, y, l, u, HX_ARRAY_SLICER_STORE)

int hx_array_vector_slicer (hx_array *x, hx_array *y,
                            int k, int loc, int dir);

#define hx_array_slice_vector(x, y, k, loc) \
  hx_array_vector_slicer(x, y, k, loc, HX_ARRAY_SLICER_SLICE)

#define hx_array_store_vector(x, y, k, loc) \
  hx_array_vector_slicer(x, y, k, loc, HX_ARRAY_SLICER_STORE)

int hx_array_matrix_slicer (hx_array *x, hx_array *y,
                            int k1, int k2, int loc,
                            int dir);

#define hx_array_slice_matrix(x, y, k1, k2, loc) \
  hx_array_matrix_slicer(x, y, k1, k2, loc, HX_ARRAY_SLICER_SLICE)

#define hx_array_store_matrix(x, y, k1, k2, loc) \
  hx_array_matrix_slicer(x, y, k1, k2, loc, HX_ARRAY_SLICER_STORE)

/* function declarations (hx-array-tile.c): */

int hx_array_tiler (hx_array *x, int k, int *nt, int *szt,
                    int dir, int incr);

#define hx_array_linearize(x, k, n, s) \
  hx_array_tiler(x, k, n, s, \
    HX_ARRAY_TILER_FORWARD, \
    HX_ARRAY_INCR_FORWARD)

#define hx_array_tileize(x, k, n, s) \
  hx_array_tiler(x, k, n, s, \
    HX_ARRAY_TILER_REVERSE, \
    HX_ARRAY_INCR_FORWARD)

int hx_array_tiling (hx_array *x, unsigned int nwords, int *nt, int *szt);

/* function declarations (hx-array-foreach.c): */

int hx_array_foreach_vector (hx_array *x, int k,
                             hx_array_foreach_cb fn, ...);

int hx_array_foreach_matrix (hx_array *x, int k1, int k2,
                             hx_array_foreach_cb fn, ...);

int hx_array_projector (hx_array *x, int k, hx_array_projector_cb fn,
                        hx_array *xp);

#endif /* __HXND_HX_ARRAY_H__ */

