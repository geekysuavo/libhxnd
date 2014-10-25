
/* ndmath: A framework for n-dimensional hypercomplex calculations for NMR.
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

/* hx_array_alloc(): allocate a hypercomplex array structure for a given
 * dimensionality. a pointer to the array structure is required by this
 * function.
 */
int hx_array_alloc (hx_array *x, int d, int k, int *sz) {
  /* declare a few required variables. */
  int i, len;

  /* check if the specified dimensionality is supported. */
  if (d < 0)
    return 0;

  /* check if the specified array dimensionality is supported. */
  if (k < 0)
    return 0;

  /* store the dimensionalities (d, k) and number of coefficients (n). */
  x->d = d;
  x->k = k;
  x->n = 1 << d;

  /* ensure that the d-dimensional shared multiplication table has been
   * initialized, and return failure if not.
   */
  if (!(x->tbl = hx_algebras_get(d)))
    return 0;

  /* allocate the array of sizes. return failure if allocation fails. */
  x->sz = (int*) calloc(x->k, sizeof(int));
  if (x->sz == NULL)
    return 0;

  /* store the values in the size array locally. */
  for (i = 0, len = x->n; i < x->k; i++) {
    /* store the size and multiply it into the total array length. */
    x->sz[i] = sz[i];
    len *= x->sz[i];
  }

  /* allocate the array of coefficients. fail if allocation fails. */
  x->x = (real*) calloc(len, sizeof(real));
  if (x->x == NULL)
    return 0;

  /* store the array total coefficient count. */
  x->len = len;

  return 1;
}

/* hx_array_free(): de-allocate a previously allocated hypercomplex array
 * structure. a pointer to the array structure is required by this function.
 */
int hx_array_free (hx_array *x) {
  /* do not attempt to free a null pointer. */
  if (x == NULL)
    return 0;

  /* check if the coefficient array is allocated. */
  if (x->x != NULL) {
    /* yes. free its associated memory. */
    free(x->x);
    x->x = NULL;
  }

  /* de-initialize the dimensionalities and coefficient count. */
  x->d = 0;
  x->n = 0;
  x->k = 0;
  x->len = 0;

  /* check if the size array is allocated. */
  if (x->sz != NULL) {
    /* yes. free its associated memory. */
    free(x->sz);
    x->sz = NULL;
  }

  /* return success. */
  return 1;
}

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
    return 0;

  /* allocate an index array. */
  arr = hx_array_index_alloc(x->k);

  /* check that the array was allocated successfully. */
  if (!arr)
    return 0;

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

/* hx_array_print(): prints a hypercomplex multidimensional array as text.
 * @x: the array to print data from.
 * @fname: the output filename.
 */
int hx_array_print (hx_array *x, const char *fname) {
  /* declare a few required variables:
   * @fh: the file handle used for writing.
   */
  int *arr, idx, i;
  FILE *fh;

  /* open the output file. */
  if (fname)
    fh = fopen(fname, "wb");
  else
    fh = stdout;

  /* check that the file was opened. */
  if (!fh)
    return 0;

  /* allocate an index array for iteration. */
  arr = hx_array_index_alloc(x->k);

  /* check that allocation was successful. */
  if (!arr)
    return 0;

  /* iterate over the points of the array. */
  idx = 0;
  do {
    /* print the indices. */
    for (i = 0; i < x->k; i++)
      fprintf(fh, "%6d ", arr[i]);

    /* print the coefficients. */
    for (i = 0; i < x->n; i++)
      fprintf(fh, "%18.8e ", x->x[i + x->n * idx]);

    /* print a newline. */
    fprintf(fh, "\n");

    /* increment the linear index. */
    idx++;
  } while (hx_array_index_inc(x->k, x->sz, &arr));

  /* free the index array. */
  free(arr);

  /* close the output file. */
  if (fname)
    fclose(fh);

  /* return success. */
  return 1;
}

/* hx_array_save(): saves a hypercomplex multidimensional array to a file.
 * @x: the array to save data from.
 * @fname: the output filename.
 */
int hx_array_save (hx_array *x, const char *fname) {
  /* declare a few required variables:
   * @wd: header that contains all array properties.
   * @fh: file handle used for writing.
   * @i: word index in the header.
   * @k: dimension index.
   */
  uint64_t wd[64];
  FILE *fh;
  int i, k;

  /* initialize the word index. */
  i = 0;

  /* open the output file. */
  if (fname)
    fh = fopen(fname, "wb");
  else
    fh = stdout;

  /* check that the file was opened. */
  if (!fh)
    return 0;

  /* zero the header bytes. */
  memset(wd, 0, 64 * sizeof(uint64_t));

  /* store a magic number into the header: 'HXNDARRY' */
  wd[i++] = (uint64_t) 0x59525241444e5848;

  /* store the array properties into the header. */
  wd[i++] = (uint64_t) x->d;
  wd[i++] = (uint64_t) x->n;
  wd[i++] = (uint64_t) x->k;
  wd[i++] = (uint64_t) x->len;
  wd[i++] = (uint64_t) sizeof(real);

  /* store the array dimension sizes into the header. */
  for (k = 0; k < x->k; k++)
    wd[i++] = (uint64_t) x->sz[k];

  /* write the header. */
  if (fwrite(wd, sizeof(uint64_t), 64, fh) != 64)
    return 0;

  /* write the array data. */
  if (fwrite(x->x, sizeof(real), x->len, fh) != x->len)
    return 0;

  /* close the output file. */
  if (fname)
    fclose(fh);

  /* return success. */
  return 1;
}

/* hx_array_load(): loads a hypercomplex multidimensional array from a file.
 * @x: the array to load data into.
 * @fname: the input filename.
 */
int hx_array_load (hx_array *x, const char *fname) {
  /* declare a few required variables:
   * @wd: header that contains all array properties.
   * @fh: file handle used for writing.
   * @i: word index in the header.
   * @k: dimension index.
   */
  uint64_t wd[64];
  FILE *fh;
  int i, k;

  /* initialize the word index. */
  i = 1;

  /* open the input file. */
  if (fname)
    fh = fopen(fname, "rb");
  else
    fh = stdin;

  /* check that the file was opened. */
  if (!fh)
    return 0;

  /* read the header from the file. */
  if (fread(wd, sizeof(uint64_t), 64, fh) != 64)
    return 0;

  /* read the array parameters from the header. */
  x->d = (int) wd[i++];
  x->n = (int) wd[i++];
  x->k = (int) wd[i++];
  x->len = (int) wd[i++];

  /* check that the data type size matches ours. if not, fail miserably.
   */
  if (wd[i++] != sizeof(real))
    return 0;

  /* allocate memory for the dimension sizes array. */
  x->sz = (int*) calloc(x->k, sizeof(int));

  /* ensure the allocation was successful. */
  if (x->sz == NULL)
    return 0;

  /* read the dimension sizes from the header. */
  for (k = 0; k < x->k; k++)
    x->sz[k] = wd[i++];

  /* allocate memory for the array data. */
  x->x = (real*) calloc(x->len, sizeof(real));

  /* ensure the allocation was successful. */
  if (!x->x)
    return 0;

  /* read the array data from the file. */
  if (fread(x->x, sizeof(real), x->len, fh) != x->len)
    return 0;

  /* ensure that the d-dimensional shared multiplication table has been
   * initialized, and return failure if not.
   */
  if (!hx_algebras_add(x->d))
    return 0;

  /* ensure that the d-dimensional shared multiplication table has been
   * initialized, and return failure if not.
   */
  if (!(x->tbl = hx_algebras_get(x->d)))
    return 0;

  /* close the input file. */
  if (fname)
    fclose(fh);

  /* and return success. */
  return 1;
}

/* hx_array_nnzdims(): find the number of nonzero-size dimensions in an array.
 * @x: a pointer to the array to use.
 */
int hx_array_nnzdims (hx_array *x) {
  /* declare a few required variables. */
  int i, nnz;

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

/* hx_array_resize(): change the dimensionality of a hypercomplex array.
 * @x: a pointer to the array to resize.
 * @d: the new algebraic dimensionality.
 * @k: the new topological dimensionality.
 * @sz: the new size array.
 */
int hx_array_resize (hx_array *x, int d, int k, int *sz) {
  /* define a required variable. */
  int *sznew, *arr, idx, idxprev, i, n, len, ok;
  int nmin, kmax;
  real *xnew;

  /* check if the specified dimensionalities are supported. */
  if (d < 0 || k < 0)
    return 0;

  /* compute the new number of coefficients. */
  n = 1 << d;

  /* compute the sizes that are smaller, before or after. */
  nmin = (n < x->n ? n : x->n);

  /* allocate a new sizes array. */
  sznew = hx_array_index_alloc(k);

  /* allocate an array to hold the iteration indices. */
  kmax = (k > x->k ? k : x->k);
  arr = hx_array_index_alloc(kmax);
  idxprev = 0;

  /* check that the size array was successfully allocated. */
  if (!sznew || !arr)
    return 0;

  /* copy the size array. */
  memcpy(sznew, sz, k * sizeof(int));

  /* compute the final size of the new coefficients array. */
  for (i = 0, len = n; i < k; i++)
    len *= sz[i];

  /* reallocate a brand new coefficients array. */
  xnew = (real*) calloc(len, sizeof(real));

  /* check that the data array was successfully allocated. */
  if (!xnew)
    return 0;

  /* loop over the set of indices. */
  for (idx = 0; idx < len / n;) {
    /* loop over the dimensions to check if we're in bounds. */
    for (i = 0, ok = 1; i < k; i++) {
      /* check if we've exceeded one of the following:
       *  - the size of the currently indexed (new) dimension.
       *  - the dimensionality of the original array.
       *  - the size of the currently indexed (original) dimension.
       */
      if (arr[i] >= sz[i] ||
          (arr[i] && (i >= x->k || arr[i] >= x->sz[i]))) {
        /* yes. break the loop. */
        ok = 0;
        break;
      }
    }

    /* check if we're on a copiable index. */
    if (ok) {
      /* pack the indices into the old array linear index. */
      hx_array_index_pack(x->k, x->sz, arr, &idxprev);

      /* loop over the coefficients. */
      for (i = 0; i < nmin; i++)
        xnew[i + n * idx] = x->x[i + x->n * idxprev];
    }

    /* increment the multidimensional indices. */
    hx_array_index_inc(k, sznew, &arr);
    idx++;
  }

  /* free the index arrays. */
  free(arr);

  /* store the new data array. */
  free(x->x);
  x->x = xnew;

  /* store the new size array. */
  free(x->sz);
  x->sz = sznew;

  /* store the new dimensionality constants. */
  x->d = d;
  x->n = n;
  x->k = k;
  x->len = len;

  /* ensure that the d-dimensional shared multiplication table has been
   * initialized, and return failure if not.
   */
  if (!(x->tbl = hx_algebras_get(d)))
    return 0;

  /* return success. */
  return 1;
}

/* hx_array_reshape(): completely reshapes an array into a new dimensionality
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

  /* compute the new data array length. */
  for (i = 0, newlen = x->n; i < k; i++)
    newlen *= sz[i];

  /* check that the new length matches the old length. */
  if (newlen != x->len)
    return 0;

  /* reallocate the sizes array. */
  x->sz = (int*) realloc(x->sz, k * sizeof(int));

  /* check that the sizes array was reallocated. */
  if (x->sz == NULL)
    return 0;

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
   */
  int ksrc, kdest;

  /* determine the source and destination array dimension indices. */
  ksrc = x->k - 1;
  kdest = x->k;

  /* check that the source dimension is tangible. */
  if (ksrc < 0)
    return 0;

  /* check that the source dimension size is divisible by the number of
   * destination dimension points.
   */
  if (x->sz[ksrc] % ndiv)
    return 0;

  /* bump the array topological dimensionality. */
  x->k++;

  /* reallocate the size array. */
  x->sz = (int*) realloc(x->sz, x->k * sizeof(int));

  /* check that the size array was reallocated successfully. */
  if (x->sz == NULL)
    return 0;

  /* set up the sizes of the source and destination array dimensions. */
  x->sz[kdest] = x->sz[ksrc] / ndiv;
  x->sz[ksrc] = ndiv;

  /* return success. */
  return 1;
}

/* hx_array_slice(): slices a portion of an array based on @lower and @upper
 * index boundaries. the operation returns its output in a new array @y,
 * which should not be allocated prior to the slice.
 * @x: pointer to the input array.
 * @y: pointer to the output array.
 * @lower: the lower index bounds.
 * @upper: the upper index bounds.
 */
int hx_array_slice (hx_array *x, hx_array *y, int *lower, int *upper) {
  /* declare a few required variables:
   */
  int i, n, ncpy, idxi, idxo, *arri, *arro, *sznew;

  /* store the number of coefficients and bytes per hypercomplex scalar. */
  n = x->n;
  ncpy = n * sizeof(real);

  /* allocate three index arrays for use during iteration. */
  arri = hx_array_index_alloc(x->k);
  arro = hx_array_index_alloc(x->k);
  sznew = hx_array_index_alloc(x->k);

  /* check that the index arrays were allocated successfully. */
  if (!arri || !arro || !sznew)
    return 0;

  /* subtract the lower bound from the upper bound. */
  hx_array_index_diff(x->k, upper, lower, &sznew);

  /* increment each element of the difference array, resulting in
   * the array of sizes of the sliced portion of the array.
   */
  for (i = 0; i < x->k; i++)
    sznew[i]++;

  /* allocate memory for the sliced outputs, but only if the output
   * array does not match the configuration of the input array.
   */
  if (hx_array_conf_cmp(x, y) && !hx_array_alloc(y, x->d, x->k, sznew))
    return 0;

  /* iterate over the larger (input) array. */
  idxi = 0;
  do {
    /* check if the current input array index is in the slice bounds. */
    if (hx_array_index_bounded(x->k, arri, lower, upper)) {
      /* yes. compute the output array indices. */
      for (i = 0; i < x->k; i++)
        arro[i] = arri[i] - lower[i];

      /* linearize the output indices and copy the coefficient memory. */
      hx_array_index_pack(x->k, sznew, arro, &idxo);
      memcpy(y->x + n * idxo, x->x + n * idxi, ncpy);
    }

    /* incremenet the input array linear index. */
    idxi++;
  } while (hx_array_index_inc(x->k, x->sz, &arri));

  /* free the allocated index arrays. */
  free(sznew);
  free(arri);
  free(arro);

  /* return success. */
  return 1;
}

/* hx_array_slice_vector(): slice a linear section from an array, starting
 * and ending at the extents of a given dimension @k. this is *much* faster
 * for slicing vectors out of nD arrays than hx_array_slice().
 * @x: pointer to the input array.
 * @y: pointer to the output array.
 * @k: array dimension to slice along.
 * @loc: off-dimension slice origin.
 */
int hx_array_slice_vector (hx_array *x, hx_array *y, int k, int *loc) {
  /* declare a few required variables:
   * @i: a loop counter used during slicing and index setup.
   * @n: size of the sliced dimension, length of the output array.
   * @ncpy: number of bytes per hypercomplex scalar value.
   * @idx: linear index of the input array.
   * @stride: stride of the input array.
   */
  int i, n, ncpy, idx, stride;

  /* check that the slice dimension in within bounds. */
  if (k < 0 || k >= x->k)
    return 0;

  /* compute the number of bytes per scalar
   * and the size of the output array.
   */
  ncpy = x->n * sizeof(real);
  n = x->sz[k];

  /* pack the off-dimension slice origin into a linear index. */
  loc[k] = 0;
  hx_array_index_pack(x->k, x->sz, loc, &idx);

  /* compute the stride for passing along the slice dimension. */
  for (i = 0, stride = 1; i < k; i++)
    stride *= x->sz[i];

  /* allocate the output array, if its configuration does not match
   * the required configuration.
   */
  if ((y->d != x->d || y->k != 1 || y->sz[0] != n) &&
      !hx_array_alloc(y, x->d, 1, &n))
    return 0;

  /* copy the scalar values into the (vector) output array. */
  for (i = 0; i < n; i++, idx += stride)
    memcpy(y->x + y->n * i, x->x + x->n * idx, ncpy);

  /* return success. */
  return 1;
}

/* hx_array_store_vector(): store a linear section from an array, essentially
 * the reverse operation of hx_array_slice_vector().
 */
int hx_array_store_vector (hx_array *x, hx_array *y, int k, int *loc) {
  /* declare a few required variables. */
  int i, n, ncpy, idx, stride;

  /* check that the slice dimension in within bounds. */
  if (k < 0 || k >= x->k)
    return 0;

  /* compute the number of bytes per scalar
   * and the size of the output array.
   */
  ncpy = x->n * sizeof(real);
  n = x->sz[k];

  /* pack the off-dimension slice origin into a linear index. */
  loc[k] = 0;
  hx_array_index_pack(x->k, x->sz, loc, &idx);

  /* compute the stride for passing along the slice dimension. */
  for (i = 0, stride = 1; i < k; i++)
    stride *= x->sz[i];

  /* check that the array configurations are correct. */
  if (y->d != x->d || y->k != 1 || y->sz[0] != n)
    return 0;

  /* copy the scalar values into the (vector) output array. */
  for (i = 0; i < n; i++, idx += stride)
    memcpy(x->x + x->n * idx, y->x + y->n * i, ncpy);

  /* return success. */
  return 1;
}

