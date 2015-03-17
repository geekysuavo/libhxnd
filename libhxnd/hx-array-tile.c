
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
 * @incr: incrementation mode, either 0 (normal) or 1 (reverse).
 */
int hx_array_tiler (hx_array *x, int k, hx_index nt, hx_index szt,
                    int dir, int incr) {
  /* declare a few required variables:
   * @pidxi: input array packed linear index.
   * @pidxo: output array packed linear index.
   * @idx: multidimensional point index array.
   * @idxt: multidimensional tile index array.
   * @xcpy: temporary tile array.
   * @res: result of inner array incrementation.
   * @rest: result of outer array incrementation.
   */
  int n, ncpy, pidxi, pidxo, res, rest;
  hx_index idx, idxt;
  hx_array xcpy;

  /* allocate the loop index arrays. */
  idx = hx_index_alloc(k);
  idxt = hx_index_alloc(k);

  /* check that all index allocations were successful. */
  if (!idx || !idxt)
    throw("failed to allocate two sets of %d indices", k);

  /* allocate a temporary array for storing tile data. */
  if (!hx_array_copy(&xcpy, x))
    throw("failed to allocate tile array");

  /* compute the coefficient count and byte count per scalar.
   */
  n = x->n;
  ncpy = n * sizeof(real);

  /* loop over the data tiles. initialize the input linear index. */
  pidxi = 0;
  do {
    /* loop over the tile points. */
    do {
      /* pack the tiled indices into a linear index. */
      hx_index_pack_tiled(k, nt, szt, idx, idxt, &pidxo);

      /* determine the copy direction. */
      switch (dir) {
        /* forward: tiles -> linear. */
        case HX_ARRAY_TILER_FORWARD:
          /* copy coefficients from the tiled array into the linear array. */
          memcpy(x->x + pidxo * n, xcpy.x + pidxi * n, ncpy);
          break;

        /* reverse: linear -> tiles. */
        case HX_ARRAY_TILER_REVERSE:
          /* copy coefficients from the linear array into the tiled array. */
          memcpy(x->x + pidxi * n, xcpy.x + pidxo * n, ncpy);
          break;
      }

      /* increment the tile linear index. */
      pidxi++;

      /* increment the point multidimensional index. */
      res = 0;
      switch (incr) {
        /* normal: lowest dimension varies fastest. */
        case HX_ARRAY_INCR_NORMAL:
          res = hx_index_incr(k, szt, idx);
          break;

        /* reverse: highest dimension varies fastest. */
        case HX_ARRAY_INCR_REVERSE:
          res = hx_index_incr_rev(k, szt, idx);
          break;
      }
    } while (res);

    /* increment the tile multidimensional index. */
    rest = 0;
    switch (incr) {
      /* normal: lowest dimension varies fastest. */
      case HX_ARRAY_INCR_NORMAL:
        rest = hx_index_incr(k, nt, idxt);
        break;

      /* reverse: highest dimension varies fastest. */
      case HX_ARRAY_INCR_REVERSE:
        rest = hx_index_incr_rev(k, nt, idxt);
        break;
    }
  } while (rest);

  /* free the allocated tile array. */
  hx_array_free(&xcpy);

  /* free the allocated multidimensional indices. */
  hx_index_free(idx);
  hx_index_free(idxt);

  /* return success. */
  return 1;
}

/* hx_array_tiling(): determine the tile counts and sizes for a hypercomplex
 * array such that the tile size is less than a certain word count.
 * @x: source array structure pointer.
 * @nwords: maximum number of data words per tile.
 * @nt: pointer to the output array of tile counts.
 * @szt: pointer to the output array of tile sizes.
 */
int hx_array_tiling (hx_array *x, unsigned int nwords,
                     hx_index nt, hx_index szt) {
  /* declare a few required variables:
   * @k: general purpose dimension loop counter.
   * @k_div: current array dimension to subdivide.
   * @n: number of current bytes per tile.
   */
  int k, k_div, n;

  /* initialize the tile counts and sizes. */
  for (k = 0; k < x->k; k++) {
    /* begin at a single tile that spans the dimension. */
    nt[k] = 1;
    szt[k] = x->sz[k];
  }

  /* initialize the current tile to be subdivided. */
  k_div = 0;

  /* loop until the tile size is below the specified limit. */
  do {
    /* do not attempt to subdivide odd tiles further. */
    while (szt[k_div] % 2 && k_div < x->k)
      k_div++;

    /* check that we found an index to subdivide. */
    if (k_div >= x->k)
      throw("failed to identify suitable tiling");

    /* divide the currently indexed tile size. */
    szt[k_div] /= 2;
    nt[k_div] *= 2;

    /* compute the new tile byte count. */
    for (k = 0, n = 1; k < x->k; k++)
      n *= szt[k];
  } while (n > nwords);

  /* return success. */
  return 1;
}

