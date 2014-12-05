
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

/* hx_array_tiler(): map scalar values in a hypercomplex array between tiled
 * array form and linear form.
 * @x: pointer to the array to linearize or tileize.
 * @k: number of tile topological dimensions.
 * @nt: array of tile counts.
 * @szt: array of tile sizes.
 * @dir: direction, either 0 (linearize) or 1 (tileize).
 */
int hx_array_tiler (hx_array *x, int k, int *nt, int *szt, int dir) {
  /* declare a few required variables:
   * @idxi: input array linear index.
   * @idxo: output array linear index.
   * @arr: multidimensional point index array.
   * @arrt: multidimensional tile index array.
   * @xcpy: temporary tile array.
   */
  int n, ncpy, idxi, idxo, *arr, *arrt;
  hx_array xcpy;

  /* allocate the loop index arrays. */
  arr = hx_array_index_alloc(k);
  arrt = hx_array_index_alloc(k);

  /* check that all index allocations were successful. */
  if (!arr || !arrt)
    throw("failed to allocate two sets of %d indices", k);

  /* allocate a temporary array for storing tile data. */
  if (!hx_array_copy(&xcpy, x))
    throw("failed to allocate tile array");

  /* compute the coefficient count and byte count per scalar.
   */
  n = x->n;
  ncpy = n * sizeof(real);

  /* loop over the data tiles. initialize the input linear index. */
  idxi = 0;
  do {
    /* loop over the tile points. */
    do {
      /* pack the tiled indices into a linear index. */
      hx_array_index_pack_tiled(k, nt, szt, arr, arrt, &idxo);

      /* determine the copy direction. */
      switch (dir) {
        /* forward: tiles -> linear. */
        case HX_ARRAY_TILER_FORWARD:
          /* copy coefficients from the tiled array into the linear array. */
          memcpy(x->x + idxo * n, xcpy.x + idxi * n, ncpy);
          break;

        /* reverse: linear -> tiles. */
        case HX_ARRAY_TILER_REVERSE:
          /* copy coefficients from the linear array into the tiled array. */
          memcpy(x->x + idxi * n, xcpy.x + idxo * n, ncpy);
          break;
      }

      /* increment the tile linear index. */
      idxi++;
    } while (hx_array_index_incr_rev(k, szt, arr));
  } while (hx_array_index_incr_rev(k, nt, arrt));

  /* free the allocated tile array. */
  hx_array_free(&xcpy);

  /* free the allocated index arrays. */
  free(arr);
  free(arrt);

  /* return success. */
  return 1;
}

