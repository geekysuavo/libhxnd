
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

/* hx_array_set_coeff(): sets a single coefficient in an array.
 * @x: pointer to the array to modify.
 * @di: coefficient index to set.
 * @value: the value of the coefficient.
 * @...: array indices to set.
 */
int hx_array_set_coeff (hx_array *x, int di, real value, ...) {
  /* declare a few required variables. */
  int i, idx, *arr;
  va_list vl;

  /* check that the coefficient index is in bounds. */
  if (di < 0 || di >= x->n)
    throw("dimension %d out of bounds [0,%d)", di, x->n);

  /* allocate an index array. */
  arr = hx_array_index_alloc(x->k);

  /* check that the array was allocated successfully. */
  if (!arr)
    throw("failed to allocate %d indices", x->k);

  /* initialize the variable arguments list. */
  va_start(vl, value);

  /* loop over the arguments list. */
  for (i = 0; i < x->k; i++)
    arr[i] = va_arg(vl, int);

  /* free the variable arguments list. */
  va_end(vl);

  /* pack the arrayed indices into a linear index. */
  hx_array_index_pack(x->k, x->sz, arr, &idx);
  idx *= x->n;
  idx += di;

  /* store the value in the coefficients array. */
  x->x[idx] = value;

  /* free the allocated index array. */
  free(arr);

  /* return success. */
  return 1;
}

/* hx_array_shuffle_block(): shuffle the values from the first and second
 * havles of a block of coefficients to yield hypercomplex data from gradient
 * enhanced data.
 * @x: the coefficient array to shuffle.
 * @buf: preallocated array the same size as @x.
 * @d: the algebraic dimensionality of the array.
 * @n: the width of each hypercomplex scalar, in reals.
 * @len: the number of hypercomplex scalars in the array.
 * @tbl: the hypercomplex multiplication table.
 *
 * NOTE: this should only be called from hx_array_complexify().
 */
int hx_array_shuffle_block (real *x, real *buf, int d, int n, int len,
                            hx_algebra tbl) {
  /* declare a few required variables:
   * @i: main scalar element loop counter.
   * @idxR: first-half array index.
   * @idxI: second-half array index.
   * @tmp, @ph: temporary scalars.
   */
  int i, idxR, idxI;
  hx_scalar ph, tmp;

  /* allocate the temporary scalars. */
  if (!hx_scalar_alloc(&ph, d) ||
      !hx_scalar_alloc(&tmp, d))
    throw("failed to allocate temporary scalars");

  /* duplicate the block into the temporary buffer. */
  memcpy(buf, x, len * n * sizeof(real));

  /* loop over the array halves. */
  for (i = 0; i < len / 2; i++) {
    /* compute the array indices. */
    idxR = n * i;
    idxI = n * (i + len / 2);

    /* perform the raw scalar data shuffle operation. */
    if (!hx_data_shuf(buf + idxR, buf + idxI, x + idxR, x + idxI,
                      ph.x, tmp.x, d, n, tbl))
      throw("failed to perform shuffle %d", i);
  }

  /* free the temporary scalars. */
  hx_scalar_free(&ph);
  hx_scalar_free(&tmp);

  /* return success. */
  return 1;
}

/* hx_array_interlace_block(): interlace the values from the first and second
 * halves of a block of coefficients into alternating values from each block.
 * @x: the coefficient array to interlace.
 * @buf: preallocated array the same size as @x.
 * @n: the number of hypercomplex scalars in the array.
 * @w: the width of each hypercomplex scalar, in reals.
 *
 * NOTE: this should only be called from hx_array_complexify().
 */
int hx_array_interlace_block (real *x, real *buf, int n, int w) {
  /* declare a few required variables:
   * @i: main scalar element loop counter.
   * @ibufR: first-half buffer (input) index.
   * @ibufI: second-half buffer (input) index.
   * @ixR: first-half output index.
   * @ixI: second-half output index.
   * @nbytes: number of bytes per scalar.
   */
  int i, ibufR, ibufI, ixR, ixI, nbytes;

  /* compute the number of bytes to copy per scalar value. */
  nbytes = w * sizeof(real);

  /* duplicate the block into the temporary buffer. */
  memcpy(buf, x, n * nbytes);

  /* loop over the array halves. */
  for (i = 0; i < n / 2; i++) {
    /* compute the source array indices. */
    ibufR = w * i;
    ibufI = w * (i + n / 2);

    /* compute the destination array indices. */
    ixR = w * (2 * i);
    ixI = w * (2 * i + 1);

    /* copy the coefficients. */
    memcpy(x + ixR, buf + ibufR, nbytes);
    memcpy(x + ixI, buf + ibufI, nbytes);
  }

  /* return success. */
  return 1;
}

/* hx_array_complexify(): convert the highest dimension (the last index) of
 * an array to the next level of algebraic dimensionality (i.e. x->d++). for
 * example, perform:
 *
 *   (d=0, k=1, sz={2*n})       -> (d=1, k=1, sz={n})
 *   (d=1, k=2, sz={n, 2*m})    -> (d=2, k=2, sz={n, m})
 *   (d=2, k=3, sz={n, m, 2*l}) -> (d=3, k=3, sz={n, m, l})
 *
 * interlacing (d = 0 -> 1):
 *   [1, u1, ...]
 *    -> [(1, u1), ...]
 *
 * interlacing (d = 1 -> 2):
 *   [(1, u1), (u2, u1u2), ...]
 *    -> [(1, u1, u2, u1u2), ...]
 *
 * interlacing (d = 2 -> 3):
 *   [(1, u1, u2, u1u2), (u3, u1u3, u2u3, u1u2u3), ...]
 *    -> [(1, u1, u2, u1u2, u3, u1u3, u2u3, u1u2u3), ...]
 *
 * it is assumed that each alternating point in the topmost array dimension
 * contains the next set of coefficients for the extra algebraic dimensions
 * created by promotion to the next greatest dimensionality.
 *
 * arguments:
 *  @x: pointer to the array to complexify, in-place.
 *  @genh: whether to apply gradient-enhanced arithmetic.
 */
int hx_array_complexify (hx_array *x, int genh) {
  /* declare a few required variables:
   * @ktop: the (topmost) array dimension to act upon.
   * @sztop: the size of the topmost array dimension.
   * @nelem: the number of scalars in the array.
   * @nblk: the number of blocks to complex-promote.
   * @szblk: the number of scalar values per block.
   * @i: the block loop counter.
   * @j: the coefficient loop counter.
   * @buf: temporary array for each block.
   */
  int ktop, sztop, nelem, nblk, szblk, i, j;
  real *buf;

  /* treat the one-dimensional (vector) array case here. in this sole
   * instance, the coefficients are already in (r,i,r,i,...) order
   */
  if (x->k == 1) {
    /* adjust the array configuration. */
    x->d++;
    x->n *= 2;
    x->sz[0] /= 2;

    /* ensure that the d-dimensional shared multiplication table has been
     * initialized, and return failure if not.
     */
    if (!(x->tbl = hx_algebras_get(x->d)))
      throw("failed to retrieve %d-algebra", x->d);

    /* return success. */
    return 1;
  }

  /* get the topmost array dimension. */
  ktop = x->k - 1;

  /* check that the topmost array dimension is valid. */
  if (ktop < 0)
    throw("topmost dimension %d is invalid", ktop);

  /* get the size of the topmost array dimension. */
  sztop = x->sz[ktop];

  /* check that the topmost array dimension is large enough. */
  if (sztop < 2)
    throw("topmost dimension %d has insufficient size %d", ktop, sztop);

  /* check that the topmost array dimension size is divisible by two. */
  if (sztop % 2)
    throw("topmost dimension %d has odd size %d", ktop, sztop);

  /* get the number of blocks and the size of each block. */
  nelem = x->len / x->n;
  nblk = sztop / 2;
  szblk = nelem / nblk;

  /* allocate a temporary swap array for interlacing. */
  buf = (real*) malloc(szblk * x->n * sizeof(real));

  /* check that allocation succeeded. */
  if (!buf)
    throw("failed to allocate interlacing buffer");

  /* loop over the blocks. */
  for (i = 0, j = 0; i < nblk; i++, j += szblk * x->n) {
    /* check if gradient-enhanced arithmetic is required. */
    if (genh) {
      /* perform gradient-enhanced arithmetic on the block. */
      if (!hx_array_shuffle_block(x->x + j, buf, x->d, x->n, szblk, x->tbl))
        throw("failed to shuffle array block %d", i);
    }

    /* interlace the block halves. */
    if (!hx_array_interlace_block(x->x + j, buf, szblk, x->n))
      throw("failed to interlace array block %d", i);
  }

  /* free the temporary swap array. */
  free(buf);

  /* store the new dimensionality and coefficient count. */
  x->d++;
  x->n *= 2;

  /* store the new topmost size. */
  x->sz[ktop] /= 2;

  /* ensure that the d-dimensional shared multiplication table has been
   * initialized, and return failure if not.
   */
  if (!(x->tbl = hx_algebras_get(x->d)))
    throw("failed to retrieve %d-algebra", x->d);

  /* return success. */
  return 1;
}

/* hx_array_real(): reduce the algebraic dimensionality of an array by
 * removing either one or all imaginary basis elements from the coefficients.
 * @x: pointer to the array to realify.
 * @d: dimension index to remove, or negative for all.
 */
int hx_array_real (hx_array *x, int d) {
  /* declare a few required variables:
   * @i: general-purpose loop counter.
   * @ord: new dimension ordering.
   */
  int i, *ord;

  /* check if all dimensions are to be realified. */
  if (d < 0) {
    /* yes. remove all imaginary coefficients. */
    return hx_array_resize(x, 0, x->k, x->sz);
  }

  /* check that the dimension index is in bounds. */
  if (d >= x->d)
    throw("dimension index %d out of bounds (-inf,%d)", x->d);

  /* allocate a dimension ordering array. */
  ord = hx_array_index_alloc(x->d);

  /* check that allocation succeeded. */
  if (!ord)
    throw("failed to allocate %d indices", x->d);

  /* build the ordering array. */
  for (i = 0; i < x->d; i++)
    ord[i] = i;

  /* store an unsortable value in the dimension to be removed. */
  ord[d] = x->d;

  /* sort the ordering array into a valid index list. */
  hx_array_index_sort(x->d, ord);

  /* reorder the array basis elements. */
  if (!hx_array_reorder_bases(x, ord))
    throw("failed to reorder basis elements");

  /* reduce the dimensionality of the array. */
  if (!hx_array_resize(x, x->d - 1, x->k, x->sz))
    throw("failed to resize array");

  /* free the ordering array. */
  free(ord);

  /* return success. */
  return 1;
}

/* hx_array_reshape(): completely reshape an array into a new dimensionality
 * @k and size set @sz.
 * @x: pointer to the array to reshape.
 * @k: the new array dimensionality.
 * @sz: the new array dimension sizes.
 */
int hx_array_reshape (hx_array *x, int k, int *sz) {
  /* declare a few required variables:
   * @i: array dimension loop counter.
   * @newlen: new array total coefficient count.
   */
  int i, newlen;

  /* check that the dimensionality is in bounds. */
  if (k < 1)
    throw("invalid dimensionality %d", k);

  /* compute the new data array length. */
  for (i = 0, newlen = x->n; i < k; i++) {
    /* check that the size is in bounds. */
    if (sz[i] < 1)
      throw("dimension size %d (#%d) out of bounds [1,inf)", sz[i], i);

    /* multiply the size into the new total length. */
    newlen *= sz[i];
  }

  /* check that the new length matches the old length. */
  if (newlen != x->len)
    throw("new size violates array length (%d != %d)", newlen, x->len);

  /* reallocate the sizes array. */
  x->sz = (int*) realloc(x->sz, k * sizeof(int));

  /* check that the sizes array was reallocated. */
  if (x->sz == NULL)
    throw("failed to reallocate %d sizes", k);

  /* store the new array dimensionality. */
  x->k = k;

  /* store the new array sizes. */
  memcpy(x->sz, sz, x->k * sizeof(int));

  /* return success. */
  return 1;
}

/* hx_array_repack(): repacks an array by taking every @ndiv points from its
 * topmost dimension and placing them into a new (added) array dimension.
 * thanks to the way the values are packed, this requires no shuffling.
 * @x: pointer to the array to repack.
 * @ndiv: the new topmost dimension size.
 */
int hx_array_repack (hx_array *x, int ndiv) {
  /* declare a few required variables:
   * @ksrc: source array dimension index.
   * @kdst: destination array dimension index.
   */
  int ksrc, kdest;

  /* determine the source and destination array dimension indices. */
  ksrc = x->k - 1;
  kdest = x->k;

  /* check that the source dimension is tangible. */
  if (ksrc < 0)
    throw("source dimension %d is invalid", ksrc);

  /* check that the target size is valid. */
  if (ndiv < 1)
    throw("target dimension size %d out of bounds [1,inf)", ndiv);

  /* check that the source dimension size is divisible by the number of
   * destination dimension points.
   */
  if (x->sz[ksrc] % ndiv)
    throw("source dimension %d is indivisible by %d", ksrc, ndiv);

  /* bump the array topological dimensionality. */
  x->k++;

  /* reallocate the size array. */
  x->sz = (int*) realloc(x->sz, x->k * sizeof(int));

  /* check that the size array was reallocated successfully. */
  if (x->sz == NULL)
    throw("failed to reallocate %d sizes", x->k);

  /* set up the sizes of the source and destination array dimensions. */
  x->sz[kdest] = x->sz[ksrc] / ndiv;
  x->sz[ksrc] = ndiv;

  /* return success. */
  return 1;
}

/* hx_array_shift_cb(): callback function for hx_array_shift().
 *
 * args:
 *  see hx_array_foreach_cb().
 *
 * varargs:
 *  @delta: shift amount, in number of coefficients.
 *  @swp: preallocated array for trace swaps.
 */
int hx_array_shift_cb (hx_array *x, hx_array *y,
                       int *arr, int idx,
                       va_list *vl) {
  /* declare required variables:
   */
  int abdelta;

  /* extract the varargs. */
  int delta = va_arg(*vl, int);
  hx_array *swp = va_arg(*vl, hx_array*);

  /* compute the absolute shift magnitude. */
  abdelta = (delta < 0 ? -delta : delta);

  /* copy first into the temporary location. */
  memcpy(swp->x, y->x, y->len * sizeof(real));

  /* act based on the shift direction. */
  if (delta < 0) {
    /* copy the shifted segments out of the temporary location. */
    memcpy(y->x, swp->x + abdelta, (y->len - abdelta) * sizeof(real));
    memcpy(y->x + (y->len - abdelta), swp->x, abdelta * sizeof(real));
  }
  else if (delta > 0) {
    /* copy the shifted segments out of the temporary location. */
    memcpy(y->x, swp->x + (y->len - abdelta), abdelta * sizeof(real));
    memcpy(y->x + abdelta, swp->x, (y->len - abdelta) * sizeof(real));
  }

  /* return success. */
  return 1;
}

/* hx_array_shift(): circularly shifts each vector along a given array
 * topological dimension, by an integer amount of points.
 * @x: pointer to the array to manipulate.
 * @k: shift topological dimension.
 * @amount: integer shift amount.
 */
int hx_array_shift (hx_array *x, int k, int amount) {
  /* declare a few required variables:
   * @swp: array for storing swapped values.
   * @delta: reduced shift amount.
   * @n: shift dimension size.
   */
  hx_array swp;
  int delta, n;

  /* check that the shift dimension index is in bounds. */
  if (k < 0 || k >= x->k)
    throw("shift index %d out of bounds [0,%d)", k, x->k);

  /* locally store the size of the shift dimension. */
  n = x->sz[k];

  /* compute the reduced shift amount. */
  delta = (amount < 0 ? -amount : amount);
  delta = delta % n;
  delta = (amount < 0 ? -delta : delta);

  /* check if the effective shift amount is zero. */
  if (delta == 0)
    return 1;

  /* multiply the shift amount by the number of coefficients per scalar. */
  delta *= x->n;

  /* allocate the temporary array. */
  if (!hx_array_alloc(&swp, x->d, 1, &n))
    throw("failed to allocate temporary shift array");

  /* perform the per-vector shift operation. */
  if (!hx_array_foreach_vector(x, k, &hx_array_shift_cb, delta, &swp))
    throw("failed to perform shift by %d", delta);

  /* free the temporary array. */
  hx_array_free(&swp);

  /* return success. */
  return 1;
}

