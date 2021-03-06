
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

/* datum_array_infill(): internally zero-fills the core array structure of an
 * NMR datum to reflect the final Nyquist grid. if no datum dimensions are
 * nonuniformly sampled, this routine has no effect.
 * @D: pointer to the datum to manipulate.
 *
 * NOTE: this routine requires that D->dims[0].sz matches D->array.sz[0].
 */
int datum_array_infill (datum *D) {
  /* declare a few required variables:
   * @i: loop counter for sampling schedule rows.
   * @j: loop counter for sampling schedule columns.
   * @pidxi: input array packed linear index.
   * @pidxo: output array packed linear index.
   * @idx: output array unpacked indices.
   * @nnus: number of nonuniform dimensions.
   * @ncpy: number of array points per sampled trace.
   * @nbytes: number of bytes per sampled trace.
   * @sznew: output array point count.
   * @tdnew: output array unpacked sizes.
   * @d: datum dimension loop counter.
   * @anew: output array structure.
   */
  int i, j, pidxi, pidxo, nnus, ncpy, nbytes, sznew;
  hx_index tdnew, idx;
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
  tdnew = hx_index_alloc(D->nd);

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
  idx = hx_index_alloc(D->nd);

  /* check that allocation was successful. */
  if (!idx)
    throw("failed to allocate %u indices", D->nd);

  /* loop over the indices in the schedule. */
  for (i = 0; i < D->n_sched; i++) {
    /* build the current unpacked index array. */
    for (j = 0, idx[0] = 0; j < D->d_sched; j++)
      idx[j + 1] = D->sched[i * D->d_sched + j];

    /* compute the input array index. */
    pidxi = i * ncpy * D->array.n;

    /* compute the output array index. */
    hx_index_pack(D->nd, tdnew, idx, &pidxo);
    pidxo *= anew.n;

    /* copy the current trace into the new array. */
    memcpy(anew.x + pidxo, D->array.x + pidxi, nbytes);
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
  hx_index_free(tdnew);
  hx_index_free(idx);

  /* return success. */
  return 1;
}

/* datum_array_refactor(): repack, infill and deinterlace the core array
 * structure of an NMR datum until it's dimensionality and complexity agree
 * with the dimension parameter values.
 * @D: pointer to the datum to manipulate.
 */
int datum_array_refactor (datum *D) {
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
    if (d == 0 && !datum_array_infill(D))
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

/* datum_array_alloc(): allocate a new array based on the currently set
 * parameters in a datum structure.
 * @D: pointer to the datum to allocate an array into.
 */
int datum_array_alloc (datum *D) {
  /* declare a few required variables. */
  unsigned int i;
  hx_index sznew;
  int d, k;

  /* check if the array has been allocated. */
  if (D->array_alloc)
    throw("array is already allocated");

  /* get the number of topological dimensions. */
  k = (int) D->nd;

  /* allocate an array for storing array dimension sizes. */
  sznew = hx_index_alloc(k);

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
  hx_index_free(sznew);

  /* return success. */
  return 1;
}

/* datum_array_read(): read and refactor the array data of a datum structure
 * that has been initialized with a correct set of parameters.
 * @D: pointer to the datum to manipulate.
 */
int datum_array_read (datum *D) {
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
  if (D->type != DATUM_TYPE_HXND && !datum_array_refactor(D))
    throw("failed to refactor array");

  /* post-process the array using format-specific routines. */
  if (!datum_type_post(D, D->type))
    throw("failed to post-process array");

  /* return success. */
  return 1;
}

/* datum_array_free(): free an allocated array structure from an NMR datum.
 * @D: pointer to the datum to manipulate.
 */
int datum_array_free (datum *D) {
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

/* datum_array_resize(): resize each array dimension (if requested) in an
 * NMR datum, but does not change the number of dimensions.
 * @D: pointer to the datum to manipulate.
 * @sz: new size array to use for resizing.
 */
int datum_array_resize (datum *D, hx_index sz) {
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

/* datum_array_slice(): slice out a portion of the array from an NMR datum.
 * @D: pointer to the datum to manipulate.
 * @lower: lower bound index array.
 * @upper: upper bound index array.
 */
int datum_array_slice (datum *D, hx_index lower, hx_index upper) {
  /* declare a few required variables:
   * @dim: datum dimension index.
   * @k: previous topological dimension count.
   * @dnew: new number of algebraic dimensions.
   * @knew: new number of topological dimensions.
   * @drm: number of removed algebraic dimensions.
   * @krm: number of removed topological dimensions.
   * @ordd: algebraic dimension reordering array.
   * @ordk: topological dimension reordering array.
   * @arrnew: destination slice array.
   */
  int dim, d, k, dnew, knew, drm, krm, dadj, kadj;
  hx_index ordd, ordk;
  hx_array arrnew;

  /* initialize the local array values, just to be safe. */
  hx_array_init(&arrnew);

  /* slice the datum array into a local array. */
  if (!hx_array_slice(&D->array, &arrnew, lower, upper))
    throw("failed to slice datum array");

  /* replace the datum array with the local array. */
  hx_array_free(&D->array);
  hx_array_copy(&D->array, &arrnew);
  hx_array_free(&arrnew);

  /* get the current dimensionalities of the datum array. */
  d = D->array.d;
  k = D->array.k;

  /* allocate arrays to specify how to reorder the datum dimensions,
   * if the need arises.
   */
  ordd = hx_index_alloc(d ? d : 1);
  ordk = hx_index_alloc(k);

  /* check that allocation succeeded. */
  if (!ordd || !ordk)
    throw("failed to allocate %d+%d indices", d, k);

  /* initialize the new dimension counts and ordering indices. */
  dnew = dadj = drm = 0;
  knew = kadj = krm = 0;

  /* compute the new topological and algebraic dimension counts. */
  for (dim = 0; dim < D->nd; dim++) {
    /* check if the current topological dimension is disappearing. */
    if (D->array.sz[D->dims[dim].k] > 1) {
      /* increment the dimension counts. */
      dnew += (D->dims[dim].d != DATUM_DIM_INVALID ? 1 : 0);
      knew++;

      /* store a sortable index in the topological ordering. */
      ordk[kadj++] = D->dims[dim].k;

      /* store an index in the algebraic ordering. */
      if (D->dims[dim].d != DATUM_DIM_INVALID)
        ordd[dadj++] = D->dims[dim].d;

      /* store the new array dimension size. */
      D->dims[dim].sz = D->array.sz[D->dims[dim].k];

      /* adjust the topological dimension index. */
      D->dims[dim].k -= krm;

      /* adjust the algebraic dimension index. */
      if (D->dims[dim].d != DATUM_DIM_INVALID)
        D->dims[dim].d -= drm;
    }
    else {
      /* store an unsortable index in the topological ordering. */
      ordk[kadj++] = k;
      krm++;

      /* store an index in the algebraic ordering. */
      if (D->dims[dim].d != DATUM_DIM_INVALID) {
        /* store an unsortable index and increment the removed count. */
        ordd[dadj++] = d;
        drm++;
      }

      /* invalidate the current topological dimension. */
      D->dims[dim].d = DATUM_DIM_INVALID;
      D->dims[dim].k = DATUM_DIM_INVALID;
      D->dims[dim].sz = 0;
    }
  }

  /* compact zero-size array dimensions out of the array. */
  if (!hx_array_compact(&D->array))
    throw("failed to compact core array");

  /* check if the dimension count changed. */
  if (knew != k) {
    /* reorder the datum dimensions. */
    if (!datum_dims_reorder(D, ordk))
      throw("failed to reorder datum dimensions");

    /* reallocate the dimension array. */
    if (!datum_dims_realloc(D, knew))
      throw("failed to reallocate dimension array");

    /* sort the ordering array into a valid index list. */
    hx_index_sort(d, ordd);

    /* reorder the core array basis elements. */
    if (!hx_array_reorder_bases(&D->array, ordd))
      throw("failed to reorder core array basis elements");

    /* reduce the dimensionality of the core array. */
    if (!hx_array_resize(&D->array, dnew, D->array.k, D->array.sz))
      throw("failed to resize core array");
  }

  /* free the allocated multidimensional indices. */
  hx_index_free(ordd);
  hx_index_free(ordk);

  /* return success. */
  return 1;
}

/* datum_array_project(): project a dimension out of an array in a datum.
 * @D: pointer to the datum to manipulate.
 * @dim: datum dimension index for projection.
 * @projector: array projection callback to use.
 */
int datum_array_project (datum *D, int dim, hx_array_projector_cb projector) {
  /* declare a few required variables:
   * @arrnew: destination projection array.
   * @d: algebraic dimension being projected.
   * @k: topological dimension being projected.
   * @ord: datum dimension reordering array.
   */
  hx_array arrnew;
  hx_index ord;
  int i, d, k;

  /* get references to the projected array dimension indices. */
  d = D->dims[dim].d;
  k = D->dims[dim].k;

  /* allocate an index array for reordering datum and array dimensions. */
  ord = hx_index_alloc(D->nd);

  /* ensure that allocation succeeded. */
  if (!ord)
    throw("failed to allocate %u indices", D->nd);

  /* apply the projector to the datum array. */
  hx_array_init(&arrnew);
  if (!hx_array_projector(&D->array, k, projector, &arrnew))
    throw("failed to apply projection operation");

  /* replace the datum array with its projection. */
  hx_array_free(&D->array);
  hx_array_copy(&D->array, &arrnew);
  hx_array_free(&arrnew);

  /* compact zero-size array dimensions out of the array. */
  if (!hx_array_compact(&D->array))
    throw("failed to compact core array");

  /* correct the topological array indices. */
  for (i = 0; i < D->nd; i++)
    D->dims[i].k -= (D->dims[i].k > k ? 1 : 0);

  /* only correct the indices if we're losing a dimension. */
  if (d != DATUM_DIM_INVALID) {
    /* correct the algebraic array indices. */
    for (i = 0; i < D->nd; i++)
      D->dims[i].d -= (D->dims[i].d > d ? 1 : 0);
  }

  /* build the datum dimension reordering array. */
  for (i = 0; i < D->nd; i++)
    ord[i] = (i == dim ? D->nd : i);

  /* reorder the datum dimensions. */
  if (!datum_dims_reorder(D, ord))
    throw("failed to reorder datum dimensions");

  /* reallocate the dimension array. */
  if (!datum_dims_realloc(D, D->nd - 1))
    throw("failed to reallocate dimension array");

  /* only reorder the complex bases if we're losing a dimension. */
  if (d != DATUM_DIM_INVALID) {
    /* build the basis reordering array. */
    for (i = 0; i < D->array.d; i++)
      ord[i] = (i == d ? D->array.d : i);

    /* sort the ordering array. */
    hx_index_sort(D->array.d, ord);

    /* reorder the core array basis elements. */
    if (!hx_array_reorder_bases(&D->array, ord))
      throw("failed to reorder core array basis elements");

    /* reduce the dimensionality of the core array. */
    if (!hx_array_resize(&D->array, D->array.d - 1, D->array.k, D->array.sz))
      throw("failed to resize core array");
  }

  /* free the allocated index array. */
  hx_index_free(ord);

  /* return success. */
  return 1;
}

