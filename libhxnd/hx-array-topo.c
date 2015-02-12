
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

/* include the n-dimensional math header. */
#include <hxnd/hx.h>

/* hx_array_nnzdims(): find the number of nonzero-size dimensions in an array.
 * @x: a pointer to the array to use.
 */
int hx_array_nnzdims (hx_array *x) {
  /* declare a few required variables. */
  int i, nnz;

  /* return zero if a null pointer was passed. */
  if (!x)
    return 0;

  /* loop over the array dimensions, counting the number of dimensions
   * that have a 'size' larger than a single element. (in other words,
   * a 1x1x5x1 array has only one nonzero-size dimension)
   */
  for (i = 0, nnz = 0; i < x->k; i++)
    if (x->sz[i] > 1) nnz++;

  /* return the number of nonzero-size dimensions. */
  return nnz;
}

/* hx_array_is_vector(): return whether an array has linear shape.
 * @x: the array to check.
 */
int hx_array_is_vector (hx_array *x) {
  /* check if the number of nonzero-size dimensions is one. */
  return (hx_array_nnzdims(x) == 1);
}

/* hx_array_is_matrix(): return whether an array has rectangular shape.
 * @x: the array to check.
 */
int hx_array_is_matrix (hx_array *x) {
  /* check if the number of nonzero-size dimensions is two. */
  return (hx_array_nnzdims(x) == 2);
}

/* hx_array_is_cube(): return whether an array has rectangular cuboid shape.
 * @x: the array to check.
 */
int hx_array_is_cube (hx_array *x) {
  /* check if the number of nonzero-size dimensions is three. */
  return (hx_array_nnzdims(x) == 3);
}

/* hx_array_compact(): removes any topological dimensions from an array that
 * have zero size, thus 'compacting' the array as much as possible.
 * @x: pointer to the array to compact.
 */
int hx_array_compact (hx_array *x) {
  /* declare a few required variables:
   * @k: original dimension index.
   * @kadj: compacted dimension index.
   * @sznew: new topological size array.
   */
  int k, kadj, *sznew;

  /* allocate a new size array. */
  sznew = hx_array_index_alloc(x->k);

  /* check that allocation succeeded. */
  if (!sznew)
    throw("failed to allocate %d indices", x->k);

  /* loop over the current dimension set. */
  for (k = 0, kadj = 0; k < x->k; k++) {
    /* check if the current dimension has nonzero size. */
    if (x->sz[k] > 1) {
      /* store the size into the new array. */
      sznew[kadj] = x->sz[k];
      kadj++;
    }
  }

  /* check if the new number of dimensions differs from the current value. */
  if (kadj != x->k) {
    /* reshape the array using the new sizes. */
    if (!hx_array_reshape(x, kadj, sznew))
      throw("failed to reshape array");
  }

  /* free the new size array. */
  free(sznew);

  /* return success. */
  return 1;
}

