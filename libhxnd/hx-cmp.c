
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

/* hx_scalar_cmp(): compare two hypercomplex scalar values for complete
 * equality.
 */
int hx_scalar_cmp (hx_scalar *a, hx_scalar *b) {
  /* declare a required variable. */
  int i;

  /* compare the algebraic dimensionality of the scalars. */
  i = hx_scalar_dims_cmp(a, b);
  if (i != HXCMP_ID)
    return i;

  /* loop over the coefficients of the scalars. */
  for (i = 0; i < a->n; i++) {
    /* compare the currently indexed coefficients. */
    if (a->x[i] < b->x[i])
      return -HXCMP_DATA;
    else if (a->x[i] > b->x[i])
      return HXCMP_DATA;
  }

  /* return equality. */
  return HXCMP_ID;
}

/* hx_scalar_dims_cmp(): compare the algebraic dimensionalities of two
 * hypercomplex scalar values for equality.
 */
int hx_scalar_dims_cmp (hx_scalar *a, hx_scalar *b) {
  /* compare the algebraic dimensionality of the scalars. */
  if (a->d < b->d)
    return -HXCMP_DIMS;
  else if (a->d > b->d)
    return HXCMP_DIMS;

  /* return equality. */
  return HXCMP_ID;
}

/* hx_array_cmp(): compare two hypercomplex arrays for complete equality.
 */
int hx_array_cmp (hx_array *a, hx_array *b) {
  /* declare a required variable. */
  int i;

  /* compare the array dimensionalities. */
  i = hx_array_dims_cmp(a, b);
  if (i != HXCMP_ID)
    return i;

  /* compare the array topologies. */
  i = hx_array_topo_cmp(a, b);
  if (i != HXCMP_ID)
    return i;

  /* loop over the coefficients array. */
  for (i = 0; i < a->len; i++) {
    /* compare the currently indexed coefficients for equality. */
    if (a->x[i] < b->x[i])
      return -HXCMP_DATA;
    else if (a->x[i] > b->x[i])
      return HXCMP_DATA;
  }

  /* return identity. */
  return HXCMP_ID;
}

/* hx_array_dims_cmp(): compare the algebraic dimensionalities of two
 * hypercomplex arrays.
 */
int hx_array_dims_cmp (hx_array *a, hx_array *b) {
  /* compare the array dimensionalities. */
  if (a->d < b->d)
    return -HXCMP_DIMS;
  else if (a->d > b->d)
    return HXCMP_DIMS;

  /* return identity. */
  return HXCMP_ID;
}

/* hx_array_topo_cmp(): compare the topologies of two hypercomplex arrays.
 */
int hx_array_topo_cmp (hx_array *a, hx_array *b) {
  /* declare a required variable. */
  int k;

  /* compare the array dimensionalities. */
  if (a->k < b->k)
    return -HXCMP_TOPO;
  else if (a->k > b->k)
    return HXCMP_TOPO;

  /* loop over the array dimensionalities. */
  for (k = 0; k < a->k; k++) {
    /* compare the size of the currently indexed array dimension. */
    if (a->sz[k] < b->sz[k])
      return -HXCMP_SIZE;
    else if (a->sz[k] > b->sz[k])
      return HXCMP_SIZE;
  }

  /* return identity. */
  return HXCMP_ID;
}

/* hx_array_conf_cmp(): compare the dimensionalities and topologies of two
 * hypercomplex arrays, collectively called the 'configurations' of the
 * arrays, for equality.
 */
int hx_array_conf_cmp (hx_array *a, hx_array *b) {
  /* declare a required variable. */
  int i;

  /* compare the array dimensionalities. */
  i = hx_array_dims_cmp(a, b);
  if (i != HXCMP_ID)
    return i;

  /* compare the array topologies. */
  i = hx_array_topo_cmp(a, b);
  if (i != HXCMP_ID)
    return i;

  /* return identity. */
  return HXCMP_ID;
}

/* hx_array_index_cmp(): compare the values in two multidimensional indices
 * of the same size @k.
 * @k: the lengths of the two index arrays.
 * @sza: the first index array.
 * @szb: the second index array.
 */
int hx_array_index_cmp (int k, int *sza, int *szb) {
  /* declare a required variable. */
  int i;

  /* return if either of the index arrays is unallocated. */
  if (!sza || !szb)
    return HXCMP_SIZE;

  /* loop over the indices. */
  for (i = 0; i < k; i++) {
    /* check if the current dimension is a mismatch. */
    if (sza[i] < szb[i])
      return -HXCMP_SIZE;
    else if (sza[i] > szb[i])
      return HXCMP_SIZE;
  }

  /* return identity. */
  return HXCMP_ID;
}

