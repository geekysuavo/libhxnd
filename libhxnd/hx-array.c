
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

/* include the byte-level data header. */
#include <hxnd/bytes.h>

/* define the number of (u64) members in the header of binary array files.
 */
#define HX_ARRAY_FWRITE_SZ_HDR  6

/* hx_array_alloc(): allocate a hypercomplex array structure for a given
 * dimensionality. a pointer to the array structure is required by this
 * function.
 */
int hx_array_alloc (hx_array *x, int d, int k, int *sz) {
  /* declare a few required variables. */
  int i, len;

  /* check if the specified dimensionality is supported. */
  if (d < 0)
    throw("invalid algebraic dimensionality %d", d);

  /* check if the specified array dimensionality is supported. */
  if (k < 1)
    throw("invalid topological dimensionality %d", k);

  /* store the dimensionalities (d, k) and number of coefficients (n). */
  x->d = d;
  x->k = k;
  x->n = 1 << d;

  /* ensure that the d-dimensional shared multiplication table has been
   * initialized, and return failure if not.
   */
  if (!(x->tbl = hx_algebras_get(d)))
    throw("failed to retrieve %d-algebra", d);

  /* allocate the array of sizes. return failure if allocation fails. */
  x->sz = hx_array_index_alloc(x->k);
  if (x->sz == NULL)
    throw("failed to allocate size array");

  /* store the values in the size array locally. */
  for (i = 0, len = x->n; i < x->k; i++) {
    /* store the size and multiply it into the total array length. */
    x->sz[i] = sz[i];
    len *= x->sz[i];
  }

  /* allocate the array of coefficients. fail if allocation fails. */
  x->x = (real*) calloc(len, sizeof(real));
  if (x->x == NULL)
    throw("failed to allocate coefficient array");

  /* store the array total coefficient count. */
  x->len = len;

  /* return success. */
  return 1;
}

/* hx_array_copy(): duplicates the contents of a hypercomplex array structure
 * into another structure.
 */
int hx_array_copy (hx_array *dst, hx_array *src) {
  /* allocate the destination array. */
  if (!hx_array_alloc(dst, src->d, src->k, src->sz))
    throw("failed to allocate destination array");

  /* copy the coefficients from the source array into the destination array.
   */
  memcpy(dst->x, src->x, src->len * sizeof(real));

  /* return success. */
  return 1;
}

/* hx_array_free(): de-allocate a previously allocated hypercomplex array
 * structure. a pointer to the array structure is required by this function.
 */
void hx_array_free (hx_array *x) {
  /* do not attempt to free a null pointer. */
  if (x == NULL)
    return;

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
    throw("failed to open '%s'", fname);

  /* allocate an index array for iteration. */
  arr = hx_array_index_alloc(x->k);

  /* check that allocation was successful. */
  if (!arr)
    throw("failed to allocate %d indices", x->k);

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
  } while (hx_array_index_incr(x->k, x->sz, arr));

  /* free the index array. */
  free(arr);

  /* close the output file. */
  if (fname)
    fclose(fh);

  /* return success. */
  return 1;
}

/* hx_array_check_magic(): checks the first "magic" bytes of a file and
 * returns whether they match the hypercomplex array magic number.
 * @fname: the input filename.
 */
int hx_array_check_magic (const char *fname) {
  /* declare a few required variables:
   * @fh: input file handle.
   * @wd: first word of input file.
   */
  uint64_t wd;
  FILE *fh;

  /* open the input file. */
  fh = fopen(fname, "rb");

  /* check that the file was opened. */
  if (!fh)
    throw("failed to open '%s'", fname);

  /* read the first word. */
  if (!fread(&wd, sizeof(uint64_t), 1, fh))
    throw("failed to read magic number");

  /* close the input file. */
  fclose(fh);

  /* check the magic word, without swapping. */
  if (wd == HX_ARRAY_MAGIC)
    return 1;

  /* swap the word and check again. */
  bytes_swap_u64(&wd);
  if (wd == HX_ARRAY_MAGIC)
    return 1;

  /* no match. */
  return 0;
}

/* hx_array_fwrite(): writes a hypercomplex multidimensional array to an
 * opened file stream.
 * @x: pointer to the source array.
 * @fh: the output file stream.
 */
int hx_array_fwrite (hx_array *x, FILE *fh) {
  /* declare a few required variables:
   * @wd: header that contains all array properties.
   * @n_wd: number of words in the header.
   * @i: word index in the header.
   * @k: dimension index.
   */
  unsigned int n_wd;
  uint64_t *wd;
  int i, k;

  /* allocate memory for the header bytes. */
  n_wd = HX_ARRAY_FWRITE_SZ_HDR + x->k;
  wd = (uint64_t*) calloc(n_wd, sizeof(uint64_t));

  /* check that the array was allocated. */
  if (!wd)
    throw("failed to allocate %d-word header", n_wd);

  /* initialize the word index. */
  i = 0;

  /* store a magic number into the header: 'HXNDARRY' */
  wd[i++] = (uint64_t) HX_ARRAY_MAGIC;

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
  if (fwrite(wd, sizeof(uint64_t), n_wd, fh) != n_wd)
    throw("failed to write %d header words", n_wd);

  /* write the array data. */
  if (fwrite(x->x, sizeof(real), x->len, fh) != x->len)
    throw("failed to write %d reals", x->len);

  /* free the header array. */
  free(wd);

  /* return success. */
  return 1;
}

/* hx_array_fread(): reads a hypercomplex multidimensional array from an
 * opened file stream.
 * @x: pointer to the destination array.
 * @fh: the input file stream.
 */
int hx_array_fread (hx_array *x, FILE *fh) {
  /* declare a few required variables:
   * @wd: header that contains all array properties.
   * @i: word index in the header.
   * @k: dimension index.
   */
  unsigned int swapping;
  uint64_t *wd1, wd0[HX_ARRAY_FWRITE_SZ_HDR];
  int i, k, n_read;

  /* read the first six words from the file. */
  n_read = fread(wd0, sizeof(uint64_t), HX_ARRAY_FWRITE_SZ_HDR, fh);
  if (n_read != HX_ARRAY_FWRITE_SZ_HDR)
    throw("failed to read initial header words");

  /* check the first word in the header. if it does not match, then
   * byte-swap the data.
   */
  if (wd0[0] != HX_ARRAY_MAGIC) {
    /* no match. swap the bytes of each word. */
    bytes_swap((uint8_t*) wd0, HX_ARRAY_FWRITE_SZ_HDR, sizeof(uint64_t));

    /* now check the magic word. */
    if (wd0[0] != HX_ARRAY_MAGIC)
      throw("invalid magic number 0x%016lx", wd0[0]);

    /* set the swapping flag. */
    swapping = 1;
  }
  else {
    /* match. no swaps needed. */
    swapping = 0;
  }

  /* initialize the word index. */
  i = 1;

  /* read the array parameters from the header. */
  x->d = (int) wd0[i++];
  x->n = (int) wd0[i++];
  x->k = (int) wd0[i++];
  x->len = (int) wd0[i++];

  /* check that the data type size matches ours. if not, fail miserably.
   */
  if (wd0[i++] != sizeof(real))
    throw("word size mismatch (%u != %u)", wd0[i - 1], sizeof(real));

  /* allocate memory for the rest of the header. */
  wd1 = (uint64_t*) calloc(x->k, sizeof(uint64_t));

  /* check that the array memory was allocated. */
  if (!wd1)
    throw("failed to allocate %d-word header", x->k);

  /* read the entire header from the file. */
  if (fread(wd1, sizeof(uint64_t), x->k, fh) != x->k)
    throw("failed to read %d header words", x->k);

  /* byte-swap, if required. */
  if (swapping)
    bytes_swap((uint8_t*) wd1, x->k, sizeof(uint64_t));

  /* allocate memory for the dimension sizes array. */
  x->sz = (int*) calloc(x->k, sizeof(int));

  /* ensure the allocation was successful. */
  if (x->sz == NULL)
    throw("failed to allocate %d sizes", x->k);

  /* read the dimension sizes from the header. */
  for (k = 0, i = 0; k < x->k; k++)
    x->sz[k] = wd1[i++];

  /* allocate memory for the array data. */
  x->x = (real*) calloc(x->len, sizeof(real));

  /* ensure the allocation was successful. */
  if (!x->x)
    throw("failed to allocate %d reals", x->len);

  /* read the array data from the file. */
  if (fread(x->x, sizeof(real), x->len, fh) != x->len)
    throw("failed to read %d reals", x->len);

  /* byte-swap, if required. */
  if (swapping)
    bytes_swap((uint8_t*) x->x, x->len, sizeof(real));

  /* ensure that the d-dimensional shared multiplication table has been
   * initialized, and return failure if not.
   */
  if (!(x->tbl = hx_algebras_get(x->d)))
    throw("failed to retrieve %d-algebra", x->d);

  /* free the header array. */
  free(wd1);

  /* and return success. */
  return 1;
}

/* hx_array_save(): saves a hypercomplex multidimensional array to a file,
 * or standard output if a NULL filename was passed.
 * @x: the array to save data from.
 * @fname: the output filename.
 */
int hx_array_save (hx_array *x, const char *fname) {
  /* declare a required variable:
   * @fh: file handle used for writing.
   */
  FILE *fh;

  /* open the output file. */
  if (fname)
    fh = fopen(fname, "wb");
  else
    fh = stdout;

  /* check that the file was opened. */
  if (!fh)
    throw("failed to open '%s'", fname);

  /* write the data to the file. */
  if (!hx_array_fwrite(x, fh))
    throw("failed to write '%s'", fname);

  /* close the output file. */
  if (fname)
    fclose(fh);

  /* return success. */
  return 1;
}

/* hx_array_load(): loads a hypercomplex multidimensional array from a file,
 * or standard input if a NULL filename is passed.
 * @x: the array to load data into.
 * @fname: the input filename.
 */
int hx_array_load (hx_array *x, const char *fname) {
  /* declare a required variable:
   * @fh: file handle used for writing.
   */
  FILE *fh;

  /* open the input file. */
  if (fname)
    fh = fopen(fname, "rb");
  else
    fh = stdin;

  /* check that the file was opened. */
  if (!fh)
    throw("failed to open '%s'", fname);

  /* read the contents of the file. */
  if (!hx_array_fread(x, fh))
    throw("failed to read '%s'", fname);

  /* close the input file. */
  if (fname)
    fclose(fh);

  /* return success. */
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

/* hx_array_interlace_block(): interlace the values from the first and second
 * halves of a block of coefficients into alternating values from each block.
 * @x: the coefficient array to interlace.
 * @buf: preallocated array the same size as @x.
 * @n: the number of hypercomplex scalars in the array.
 * @w: the width of each hypercomplex scalar, in reals.
 */
int hx_array_interlace_block (real *x, real *buf, int n, int w) {
  /* declare a few required variables:
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
    /* FIXME: handle gradient-enhanced arithmetic on a per-block basis. */

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

  /* return success. */
  return 1;
}

/* hx_array_resize(): change the configuration of a hypercomplex array.
 * @x: a pointer to the array to resize.
 * @d: the new algebraic dimensionality.
 * @k: the new topological dimensionality.
 * @sz: the new size array.
 */
int hx_array_resize (hx_array *x, int d, int k, int *sz) {
  /* define a required variable. */
  int *arr, idx, idxprev, i, n, len, ok;
  int nmin, kmax;
  real *xnew;

  /* check if the specified dimensionalities are supported. */
  if (d < 0 || k < 1)
    throw("dimensions (%d, %d) are invalid", d, k);

  /* compute the new number of coefficients. */
  n = 1 << d;

  /* compute the sizes that are smaller, before or after. */
  nmin = (n < x->n ? n : x->n);

  /* allocate an array to hold the iteration indices. */
  kmax = (k > x->k ? k : x->k);
  arr = hx_array_index_alloc(kmax);
  idxprev = 0;

  /* check that the index array was successfully allocated. */
  if (!arr)
    throw("failed to allocate %d indices", kmax);

  /* compute the final size of the new coefficients array. */
  for (i = 0, len = n; i < k; i++)
    len *= sz[i];

  /* reallocate a brand new coefficients array. */
  xnew = (real*) calloc(len, sizeof(real));

  /* check that the data array was successfully allocated. */
  if (!xnew)
    throw("failed to allocate %d reals", len);

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
    hx_array_index_incr(k, sz, arr);
    idx++;
  }

  /* free the index arrays. */
  free(arr);

  /* store the new data array. */
  free(x->x);
  x->x = xnew;

  /* resize the size array, if the new size is different. */
  if (k != x->k)
    x->sz = (int*) realloc(x->sz, k * sizeof(int));

  /* store the new size array. */
  if (sz != x->sz)
    memcpy(x->sz, sz, k * sizeof(int));

  /* store the new dimensionality constants. */
  x->d = d;
  x->n = n;
  x->k = k;
  x->len = len;

  /* ensure that the d-dimensional shared multiplication table has been
   * initialized, and return failure if not.
   */
  if (!(x->tbl = hx_algebras_get(d)))
    throw("failed to retrieve %d-algebra", d);

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

/* hx_array_compact(): removes any topological dimensions from an array that
 * have zero size, thus 'compacting' the array as much as possible.
 * @x: pointer to the array to compact.
 */
int hx_array_compact (hx_array *x) {
  /* declare a few required variables:
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

/* hx_array_slicer(): slices a portion of an array based on @lower and @upper
 * index boundaries. the operation returns its output in a new array @y,
 * which should not be allocated prior to the slice.
 * @x: pointer to the input array.
 * @y: pointer to the output array.
 * @lower: the lower index bounds.
 * @upper: the upper index bounds.
 * @dir: either HX_ARRAY_SLICER_SLICE or HX_ARRAY_SLICER_STORE.
 */
int hx_array_slicer (hx_array *x, hx_array *y,
                     int *lower, int *upper,
                     int dir) {
  /* declare a few required variables:
   * @i: general-purpose loop counter.
   * @xycmp: set high if @y needs allocation.
   * @n: number of coefficients per scalar.
   * @ncpy: number of bytes per scalar.
   * @idxi: input array linear index.
   * @idxo: output array linear index.
   * @arri: input array index set.
   * @arro: output array index set.
   * @sznew: output array sizes.
   */
  int i, xycmp, n, ncpy, idxi, idxo, *arri, *arro, *sznew;

  /* store the number of coefficients and bytes per hypercomplex scalar. */
  n = x->n;
  ncpy = n * sizeof(real);

  /* allocate three index arrays for use during iteration. */
  arri = hx_array_index_alloc(x->k);
  arro = hx_array_index_alloc(x->k);
  sznew = hx_array_index_alloc(x->k);

  /* check that the index arrays were allocated successfully. */
  if (!arri || !arro || !sznew)
    throw("failed to allocate %d indices", x->k);

  /* subtract the lower bound from the upper bound. */
  hx_array_index_diff(x->k, upper, lower, sznew);

  /* increment each element of the difference array, resulting in
   * the array of sizes of the sliced portion of the array.
   */
  for (i = 0; i < x->k; i++)
    sznew[i]++;

  /* check that the output array has the same dimensionalities as
   * the input array.
   */
  xycmp = ((x->d != y->d) || (x->k != y->k));

  /* if the above comparisons both evaluated to false, then the following
   * block is safe to execute.
   */
  if (!xycmp) {
    /* check that the output array has the desired size. */
    for (i = 0; i < x->k; i++)
      xycmp = (xycmp || (y->sz[i] != sznew[i]));
  }

  /* allocate memory for the sliced outputs, but only if the output
   * array does pass the tests performed above that yield @xycmp.
   */
  if (xycmp && !hx_array_alloc(y, x->d, x->k, sznew))
    throw("failed to allocate slice destination array");

  /* iterate over the larger (input) array. */
  idxi = 0;
  do {
    /* check if the current input array index is in the slice bounds. */
    if (hx_array_index_bounded(x->k, arri, lower, upper)) {
      /* yes. compute the output array indices. */
      for (i = 0; i < x->k; i++)
        arro[i] = arri[i] - lower[i];

      /* linearize the output indices. */
      hx_array_index_pack(x->k, sznew, arro, &idxo);

      /* copy the coefficient memory. */
      switch (dir) {
        /* slice: x ==> y */
        case HX_ARRAY_SLICER_SLICE:
          memcpy(y->x + n * idxo, x->x + n * idxi, ncpy);
          break;

        /* store: x <== y */
        case HX_ARRAY_SLICER_STORE:
          memcpy(x->x + n * idxi, y->x + n * idxo, ncpy);
          break;

        /* other: no-op. */
        default:
          break;
      }
    }

    /* incremenet the input array linear index. */
    idxi++;
  } while (hx_array_index_incr(x->k, x->sz, arri));

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
int hx_array_slice_vector (hx_array *x, hx_array *y, int k, int loc) {
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
    throw("slice dimension %d out of bounds [0,%d)", k, x->k);

  /* compute the number of bytes per scalar
   * and the size of the output array.
   */
  ncpy = x->n * sizeof(real);
  n = x->sz[k];

  /* compute the stride for passing along the slice dimension. */
  for (i = 0, stride = 1; i < k; i++)
    stride *= x->sz[i];

  /* allocate the output array, if its configuration does not match
   * the required configuration.
   */
  if ((y->d != x->d || y->k != 1 || y->sz[0] != n) &&
      !hx_array_alloc(y, x->d, 1, &n))
    throw("failed to allocate slice destination array");

  /* copy the scalar values into the (vector) output array. */
  for (i = 0, idx = loc; i < n; i++, idx += stride)
    memcpy(y->x + y->n * i, x->x + x->n * idx, ncpy);

  /* return success. */
  return 1;
}

/* hx_array_store_vector(): store a linear section from an array, essentially
 * the reverse operation of hx_array_slice_vector().
 */
int hx_array_store_vector (hx_array *x, hx_array *y, int k, int loc) {
  /* declare a few required variables. */
  int i, n, ncpy, idx, stride;

  /* check that the slice dimension in within bounds. */
  if (k < 0 || k >= x->k)
    throw("slice dimension %d out of bounds [0,%d)", k, x->k);

  /* compute the number of bytes per scalar
   * and the size of the output array.
   */
  ncpy = x->n * sizeof(real);
  n = x->sz[k];

  /* compute the stride for passing along the slice dimension. */
  for (i = 0, stride = 1; i < k; i++)
    stride *= x->sz[i];

  /* check that the array configurations are correct. */
  if (y->d != x->d || y->k != 1 || y->sz[0] != n)
    throw("source-destination array configuration mismatch");

  /* copy the scalar values into the (vector) output array. */
  for (i = 0, idx = loc; i < n; i++, idx += stride)
    memcpy(x->x + x->n * idx, y->x + y->n * i, ncpy);

  /* return success. */
  return 1;
}

/* hx_array_vector_op(): perform an operation on each vector of a hypercomplex
 * array @x using a standardized callback-based scheme.
 */
int hx_array_vector_op (hx_array *x, int k, hx_array_vector_cb fn, ...) {
  /* declare a few required variables:
   * @arr: the index array for the current position in @x.
   * @idx: the linear index for the current position in @x.
   * @y: hypercomplex array holding the currently sliced vector values.
   * @vl: the variable argument list passed to each callback invocation.
   */
  int *arr, idx, szk, slice;
  hx_array y;
  va_list vl;

  /* check that the dimension index is in bounds. */
  if (k < 0 || k >= x->k)
    throw("dimension index %d is out of bounds [0,%d)", k, x->k);

  /* retrieve the size of the dimension under operation. */
  szk = x->sz[k];

  /* allocate an index array. */
  arr = hx_array_index_alloc(x->k);

  /* check that allocation was successful. */
  if (!arr)
    throw("failed to allocate %d indices", x->k);

  /* allocate a temporary array to store each sliced vector. */
  if (!hx_array_alloc(&y, x->d, 1, &szk))
    throw("failed to allocate slice (%d, 1)-array", x->d);

  /* initialize the linear indices. */
  idx = slice = 0;

  /* iterate over the elements of the array. */
  do {
    /* pack the index array into a linear index. */
    hx_array_index_pack(x->k, x->sz, arr, &idx);

    /* slice the currently indexed vector from the array. */
    if (!hx_array_slice_vector(x, &y, k, idx))
      throw("failed to slice vector %d", slice);

    /* initialize the variable arguments list. */
    va_start(vl, fn);

    /* execute the callback function. */
    if (!fn(x, &y, arr, idx, &vl))
      throw("failed to execute callback %d", slice);

    /* free the variable arguments list. */
    va_end(vl);

    /* store the modified sliced vector back into the array. */
    if (!hx_array_store_vector(x, &y, k, idx))
      throw("failed to store vector %d", slice);

    /* increment the slice index. */
    slice++;
  }
  while (hx_array_index_skip(x->k, x->sz, arr, k));

  /* free the temporary array. */
  hx_array_free(&y);

  /* free the index array. */
  free(arr);

  /* return success. */
  return 1;
}

/* hx_array_shift_cb(): callback function for hx_array_shift().
 *
 * args:
 *  see hx_array_vector_cb().
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
  if (!hx_array_vector_op(x, k, &hx_array_shift_cb, delta, &swp))
    throw("failed to perform shift by %d", delta);

  /* free the temporary array. */
  hx_array_free(&swp);

  /* return success. */
  return 1;
}

