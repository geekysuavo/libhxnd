
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

/* include the nmr data header. */
#include <hxnd/nmr-datum.h>

/* declare a buffer size for reading bruker/varian schedule files.
 */
#define N_BUF  256

/* datum_infill_array(): internally zero-fills the core array structure of an
 * NMR datum to reflect the final Nyquist grid. if no datum dimensions are
 * nonuniformly sampled, this routine has no effect.
 * @D: pointer to the datum to manipulate.
 *
 * NOTE: this routine requires that D->dims[0].sz matches D->array.sz[0].
 */
int datum_infill_array (datum *D) {
  /* declare a few required variables:
   * @i: loop counter for sampling schedule rows.
   * @j: loop counter for sampling schedule columns.
   * @idxi: input array linear index.
   * @idxo: output array linear index.
   * @arr: output array unpacked indices.
   * @nnus: number of nonuniform dimensions.
   * @ncpy: number of array points per sampled trace.
   * @nbytes: number of bytes per sampled trace.
   * @sznew: output array point count.
   * @tdnew: output array unpacked sizes.
   * @d: datum dimension loop counter.
   * @anew: output array structure.
   */
  int i, j, idxi, idxo, nnus, ncpy, nbytes, sznew, *tdnew, *arr;
  unsigned int d;
  hx_array anew;

  /* determine the number of nonuniform dimensions. */
  for (d = 0, nnus = 0; d < D->nd; d++)
    nnus += (D->dims[d].nus ? 1 : 0);

  /* return successfully if all dimensions are uniform. */
  if (nnus == 0)
    return 1;

  /* check that the datum contains a schedule array. */
  if (D->sched == NULL || D->d_sched < 1 || D->n_sched < 1)
    throw("datum contains no schedule array");

  /* allocate a new size array to build the new array. */
  tdnew = hx_array_index_alloc(D->nd);

  /* check that the size array was allocated. */
  if (!tdnew)
    throw("failed to allocate %u indices", D->nd);

  /* initialize the calculated values. */
  sznew = ncpy = D->dims[0].sz;
  tdnew[0] = D->dims[0].td;

  /* loop over the datum dimensions to build the new size array. */
  for (d = 1; d < D->nd; d++) {
    /* compute the new time-domain point count. */
    tdnew[d] = (D->dims[d].nus ? D->dims[d].tdunif : D->dims[d].td);

    /* compute the number of points to copy per trace. */
    ncpy *= (D->dims[d].cx ? 2 : 1);

    /* compute the array length. */
    sznew *= tdnew[d];
  }

  /* compute the number of bytes to copy per trace. */
  nbytes = ncpy * D->array.n * sizeof(real);

  /* allocate a new array to store the infilled values. */
  if (!hx_array_alloc(&anew, D->array.d, D->array.k, &sznew))
    throw("failed to allocate new (%d, %d)-array", D->array.d, D->array.k);

  /* allocate a new index array for traversal. */
  arr = hx_array_index_alloc(D->nd);

  /* check that allocation was successful. */
  if (!arr)
    throw("failed to allocate %u indices", D->nd);

  /* loop over the indices in the schedule. */
  for (i = 0; i < D->n_sched; i++) {
    /* build the current unpacked index array. */
    for (j = 0, arr[0] = 0; j < D->d_sched; j++)
      arr[j + 1] = D->sched[i * D->d_sched + j];

    /* compute the input array index. */
    idxi = i * ncpy * D->array.n;

    /* compute the output array index. */
    hx_array_index_pack(D->nd, tdnew, arr, &idxo);
    idxo *= anew.n;

    /* copy the current trace into the new array. */
    memcpy(anew.x + idxo, D->array.x + idxi, nbytes);
  }

  /* replaced the datum array with the local infilled array. */
  hx_array_free(&D->array);
  hx_array_copy(&D->array, &anew);
  hx_array_free(&anew);

  /* store the new sizes in the datum dimensions. */
  for (d = 1; d < D->nd; d++) {
    /* store the time-domain size. */
    D->dims[d].sz = D->dims[d].td = tdnew[d];

    /* check if the dimension is complex. */
    if (D->dims[d].cx)
      D->dims[d].sz /= 2;
  }

  /* free the allocated index arrays. */
  free(tdnew);
  free(arr);

  /* return success. */
  return 1;
}

/* datum_refactor_array(): repacks, infills and deinterlaces the core array
 * structure of an NMR datum until it's dimensionality and complexity agree
 * with the dimension parameter values.
 * @D: pointer to the datum to manipulate.
 */
int datum_refactor_array (datum *D) {
  /* declare a few required variables:
   * @d: dimension loop counter.
   */
  unsigned int d;

  /* check that the array has been allocated. */
  if (!D->array_alloc)
    throw("array is unallocated");

  /* loop over the acquisition dimensions to refactor the nD array. */
  for (d = 0; d < D->nd; d++) {
    /* repack indirect dimensions in the array. */
    if (d > 0 && !hx_array_repack(&D->array, D->dims[d - 1].sz))
      throw("failed to repack array dimension %d", d);

    /* check if the current dimension is complex. */
    if (D->dims[d].cx) {
      /* complexify this dimension. */
      if (!hx_array_complexify(&D->array, D->dims[d].genh))
        throw("failed to complexify dimension %d", d);

      /* set the algebraic dimension index. */
      D->dims[d].d = D->array.d - 1;
    }
    else {
      /* set /no/ algebraic dimension index. */
      D->dims[d].d = DATUM_DIM_INVALID;
    }

    /* infill nonuniformly sampled indirect dimensions. */
    if (d == 0 && !datum_infill_array(D))
      throw("failed to infill nonuniformly sampled dimensions");

    /* set the topological dimension index. */
    D->dims[d].k = (int) d;
  }

  /* handle post-refactor array tweaking operations. */
  for (d = 0; d < D->nd; d++) {
    /* check if sign alternation is needed. */
    if (D->dims[d].alt) {
      /* apply sign alternation. */
      if (!hx_array_alternate_sign(&D->array, D->dims[d].k))
        throw("failed to sign-alternate dimension %d", D->dims[d].k);
    }

    /* check if imaginary negation is needed. */
    if (D->dims[d].neg) {
      /* apply imaginary negation. */
      if (!hx_array_negate_basis(&D->array, D->dims[d].d))
        throw("failed to negate imaginary dimension %d", D->dims[d].d);
    }
  }

  /* return success. */
  return 1;
}

/* datum_read_sched(): reads a bruker/varian schedule file into the schedule
 * array of a datum structure.
 * @D: pointer to the datum structure to manipulate.
 * @fname: the input filename.
 */
int datum_read_sched (datum *D, const char *fname) {
  /* declare required variables:
   * @buf: character buffer for reading lines.
   * @tokv: array of token strings.
   * @tokc: number of strings in @tokc.
   * @d: number of schedule columns.
   * @n: number of schedule rows.
   * @i: general purpose loop index.
   * @sched: output schedule array.
   * @fh: input file handle.
   */
  char buf[N_BUF], **tokv;
  unsigned int toki, tokc;
  int d, n, i, *sched;
  FILE *fh;

  /* initialize the results. */
  sched = NULL;
  d = n = 0;

  /* open the input file. */
  fh = fopen(fname, "rb");

  /* check that the input file was opened. */
  if (!fh)
    throw("failed to open '%s'", fname);

  /* initialize the array index. */
  i = 0;

  /* loop until we've read the entire file. */
  while (!feof(fh)) {
    /* read a new line from the file. */
    if (fgets(buf, N_BUF, fh)) {
      /* split the read line into tokens. */
      tokv = strsplit(buf, " ", &tokc);

      /* check that the split was successful. */
      if (!tokv || tokc < 1)
        throw("failed to split string '%s'", buf);

      /* trim the token strings. */
      strvtrim(tokv, tokc);

      /* store/check the number of array elements. */
      if (i == 0)
        d = tokc;
      else if (tokc != d)
        throw("unexpected token count %d", tokc);

      /* increment the row count. */
      n++;

      /* reallocate the schedule array. */
      sched = (int*) realloc(sched, d * n * sizeof(int));

      /* check that the reallocation succeeded. */
      if (!sched)
        throw("failed to reallocate schedule array");

      /* loop over the tokens. */
      for (toki = 0; toki < tokc; toki++, i++)
        sched[i] = atoi(tokv[toki]);

      /* free the string array. */
      strvfree(tokv, tokc);
    }
  }

  /* close the input file. */
  fclose(fh);

  /* store the identified array parameters. */
  D->d_sched = d;
  D->n_sched = n;

  /* store the constructed schedule array. */
  D->sched = sched;

  /* return success. */
  return 1;
}

/* datum_alloc_array(): allocates a new array based on the currently set
 * parameters in a datum structure.
 * @D: pointer to the datum to allocate an array into.
 */
int datum_alloc_array (datum *D) {
  /* declare a few required variables. */
  int d, k, *sznew;
  unsigned int i;

  /* check if the array has been allocated. */
  if (D->array_alloc)
    throw("array is already allocated");

  /* get the number of topological dimensions. */
  k = (int) D->nd;

  /* allocate an array for storing array dimension sizes. */
  sznew = hx_array_index_alloc(k);

  /* check that allocation was successful. */
  if (!sznew)
    throw("failed to allocate %d indices", k);

  /* compute the number of algebraic dimensions. */
  for (i = 0, d = 0; i < D->nd; i++) {
    /* store the topological dimension index. */
    D->dims[i].k = (int) i;

    /* store the topological dimension size. */
    sznew[i] = (D->dims[i].sz ? D->dims[i].sz : 1);

    /* check if the dimension is complex. */
    if (D->dims[i].cx) {
      /* store the algebraic dimension index. */
      D->dims[i].d = d++;
    }
    else {
      /* store an invalid dimension index. */
      D->dims[i].d = DATUM_DIM_INVALID;
    }
  }

  /* allocate the core array. */
  if (!hx_array_alloc(&D->array, d, k, sznew))
    throw("failed to allocate core array");

  /* indicate that the array has been allocated. */
  D->array_alloc = 1;

  /* free the size array. */
  free(sznew);

  /* return success. */
  return 1;
}

/* datum_read_array(): read and refactor the array data of a datum structure
 * that has been initialized with a correct set of parameters.
 * @D: pointer to the datum to manipulate.
 */
int datum_read_array (datum *D) {
  /* check if the array has been allocated. */
  if (D->array_alloc)
    return 1;

  /* check that the filename is non-null and non-empty. */
  if (D->fname == NULL || strlen(D->fname) == 0)
    throw("filename is invalid");

  /* load the array data using format-specific routines. */
  if (!datum_type_array(D, D->type))
    throw("failed to read array data");

  /* indicate that the array has been allocated. */
  D->array_alloc = 1;

  /* refactor the core array, but only if it's *not* already in
   * the native hxnd format.
   */
  if (D->type != DATUM_TYPE_HXND && !datum_refactor_array(D))
    throw("failed to refactor array");

  /* post-process the array using format-specific routines. */
  if (!datum_type_post(D, D->type))
    throw("failed to post-process array");

  /* return success. */
  return 1;
}

/* datum_free_array(): frees an allocated array structure from an NMR datum.
 * @D: pointer to the datum to manipulate.
 */
int datum_free_array (datum *D) {
  /* check if the array is allocated. */
  if (!D->array_alloc)
    return 1;

  /* free the data array. */
  hx_array_free(&D->array);

  /* indicate that the array has been de-allocated. */
  D->array_alloc = 0;

  /* return success. */
  return 1;
}

/* datum_resize_array(): resizes each array dimension (if requested) in an
 * NMR datum, but does not change the number of dimensions.
 * @D: pointer to the datum to manipulate.
 * @sz: new size array to use for resizing.
 */
int datum_resize_array (datum *D, int *sz) {
  /* declare a few required variables:
   * @d: dimension loop counter.
   */
  unsigned int d;

  /* check that the array has been allocated. */
  if (!D->array_alloc)
    throw("array is unallocated");

  /* check that each new size is greater than one. */
  for (d = 0; d < D->nd; d++) {
    /* check that the current size is in bounds. */
    if (sz[d] < 2)
      throw("invalid size %d along dimension %u", sz[d], d);
  }

  /* attempt to resize the array content. */
  if (!hx_array_resize(&D->array, D->array.d, D->array.k, sz))
    throw("failed to resize core array");

  /* loop over the acquisition dimensions to store the new sizes. */
  for (d = 0; d < D->nd; d++)
    D->dims[d].sz = sz[d];

  /* return success. */
  return 1;
}

/* datum_slice_array(): slices out a portion of the array from an NMR datum.
 * @D: pointer to the datum to manipulate.
 * @lower: lower bound index array.
 * @upper: upper bound index array.
 */
int datum_slice_array (datum *D, int *lower, int *upper) {
  /* declare a few required variables:
   * @arrnew: destination slice array.
   * @nd: previous dimension count.
   * @ndnew: new number of dimensions.
   * @ord: dimension reordering array.
   */
  int nd, ndnew, d, dadj, *ord;
  hx_array arrnew;

  /* initialize the local array values, just to be safe. */
  arrnew.d = arrnew.k = 0;
  arrnew.sz = NULL;

  /* slice the datum array into a local array. */
  if (!hx_array_slice(&D->array, &arrnew, lower, upper))
    throw("failed to slice datum array");

  /* replace the datum array with the local array. */
  hx_array_free(&D->array);
  hx_array_copy(&D->array, &arrnew);
  hx_array_free(&arrnew);

  /* allocate an array to specify how to reorder the datum dimensions,
   * if the need arises.
   */
  nd = D->nd;
  ord = hx_array_index_alloc(nd);

  /* check that allocation succeeded. */
  if (!ord)
    throw("failed to allocate %u indices", nd);

  /* store the new array sizes in the datum dimensions. */
  for (d = 0, dadj = 0, ndnew = 0; d < nd; d++) {
    /* check if the current dimension has nonzero size. */
    if (D->array.sz[d] > 1) {
      /* store a sortable value in the ordering array. */
      ord[d] = dadj;
      dadj++;

      /* increment the number of nonzero-size dimensions. */
      ndnew++;
    }
    else {
      /* store an unsortable value in the ordering array. */
      ord[d] = (int) nd;
    }

    /* store the new array dimension size. */
    D->dims[d].sz = D->array.sz[d];
  }

  /* compact zero-size array dimensions out of the array. */
  if (!hx_array_compact(&D->array))
    throw("failed to compact core array");

  /* check if the dimension count changed. */
  if (ndnew != nd) {
    /* reorder the datum dimensions. */
    if (!datum_dims_reorder(D, ord))
      throw("failed to reorder datum dimensions");

    /* reallocate the dimension array. */
    if (!datum_dims_realloc(D, ndnew))
      throw("failed to reallocate dimension array");

    /* sort the ordering array into a valid index list. */
    hx_array_index_sort(nd, ord);

    /* reorder the core array basis elements. */
    if (!hx_array_reorder_bases(&D->array, ord))
      throw("failed to reorder core array basis elements");

    /* reduce the dimensionality of the core array. */
    if (!hx_array_resize(&D->array, ndnew, D->array.k, D->array.sz))
      throw("failed to resize core array");
  }

  /* free the allocated index array. */
  free(ord);

  /* return success. */
  return 1;
}

