
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

/* datum_dim_desc: structure that defines a mapping between dimension
 * parameter names, their locations within the datum_dim structure,
 * and their types in the structure.
 */
typedef struct {
  /* @name: parameter name string.
   * @type: parameter type character.
   * @off: struct member offset.
   */
  const char *name;
  char type;
  size_t off;
}
datum_dim_desc;

/* datum_dim_parms: array of dimension description structures that indicate
 * the names, types and offets of each datum dimension structure member .
 */
static const datum_dim_desc datum_dim_parms[] = {
  /* size parameters. */
  { "sz",        'u', offsetof(datum_dim, sz) },
  { "td",        'u', offsetof(datum_dim, td) },
  { "tdunif",    'u', offsetof(datum_dim, tdunif) },

  /* status flags. */
  { "complex",   'b', offsetof(datum_dim, cx) },
  { "nus",       'b', offsetof(datum_dim, nus) },
  { "ft",        'b', offsetof(datum_dim, ft) },
  { "alternate", 'b', offsetof(datum_dim, alt) },
  { "negate",    'b', offsetof(datum_dim, neg) },
  { "gradient",  'b', offsetof(datum_dim, genh) },

  /* spectral parameters. */
  { "carrier",   'f', offsetof(datum_dim, carrier) },
  { "width",     'f', offsetof(datum_dim, width) },
  { "offset",    'f', offsetof(datum_dim, offset) },

  /* nucleus string. */
  { "name",      's', offsetof(datum_dim, nuc) },

  /* null-terminator. */
  { NULL, 'c', 0 }
};

/* datum_dims_getparm(): retrieve a dimension parameter from a datum
 * structure, by name and dimension index.
 * @D: pointer to the datum struct.
 * @name: parameter name to get.
 * @d: dimension of parameter.
 * @parm: pointer to the result.
 */
int datum_dims_getparm (datum *D, const char *name,
                        unsigned int d,
                        void *parm) {
  /* declare pointer locations for each data type. */
  unsigned int *uptr;
  char *sptr;
  real *fptr;
  int *dptr;

  /* declare required variables. */
  unsigned int i;
  uint8_t *base;

  /* check that the dimension index is within bounds. */
  if (d >= D->nd)
    throw("dimension index %u out of bound %u", d, D->nd);

  /* compute the base address of the dimension structure of interest. */
  base = (uint8_t*) &D->dims[d];

  /* loop over the dimension descriptors. */
  for (i = 0; datum_dim_parms[i].name; i++) {
    /* check if the current descriptor matches the requested name. */
    if (strcmp(datum_dim_parms[i].name, name) == 0) {
      /* act based on the parameter type. */
      switch (datum_dim_parms[i].type) {
        /* signed int */
        case 'd':
          dptr = (int*) (base + datum_dim_parms[i].off);
          *((int*) parm) = *dptr;
          return 1;

        /* unsigned int (also boolean) */
        case 'u':
        case 'b':
          uptr = (unsigned int*) (base + datum_dim_parms[i].off);
          *((unsigned int*) parm) = *uptr;
          return 1;

        /* real */
        case 'f':
          fptr = (real*) (base + datum_dim_parms[i].off);
          *((real*) parm) = *fptr;
          return 1;

        /* string */
        case 's':
          sptr = (char*) (base + datum_dim_parms[i].off);
          memcpy(parm, sptr, 8);
          return 1;

        /* unknown */
        default:
          throw("invalid parameter type '%c'", datum_dim_parms[i].type);
      }
    }
  }

  /* no match was found. throw an exception. */
  throw("invalid dimension parameter name '%s'", name);
}

/* datum_dims_setparm(): store a dimension parameter into a datum
 * structure, by name and dimension index.
 * @D: pointer to the datum struct.
 * @name: parameter name to set.
 * @d: dimension of parameter.
 * @parm: the new value, as a string.
 */
int datum_dims_setparm (datum *D, const char *name,
                        unsigned int d,
                        const char *parm) {
  /* declare pointer locations for each data type. */
  unsigned int *uptr;
  char *sptr;
  real *fptr;
  int *dptr;

  /* declare required variables. */
  unsigned int i;
  uint8_t *base;

  /* check that the dimension index is within bounds. */
  if (d >= D->nd)
    throw("dimension index %u out of bound %u", d, D->nd);

  /* compute the base address of the dimension structure of interest. */
  base = (uint8_t*) &D->dims[d];

  /* loop over the dimension descriptors. */
  for (i = 0; datum_dim_parms[i].name; i++) {
    /* check if the current descriptor matches the requested name. */
    if (strcmp(datum_dim_parms[i].name, name) == 0) {
      /* act based on the parameter type. */
      switch (datum_dim_parms[i].type) {
        /* signed int */
        case 'd':
          dptr = (int*) (base + datum_dim_parms[i].off);
          *dptr = atol(parm);
          return 1;

        /* unsigned int */
        case 'u':
          uptr = (unsigned int*) (base + datum_dim_parms[i].off);
          *uptr = (unsigned int) atol(parm);
          return 1;

        /* boolean (stored as unsigned int) */
        case 'b':
          uptr = (unsigned int*) (base + datum_dim_parms[i].off);
          *uptr = (unsigned int) strbool(parm);
          return 1;

        /* real */
        case 'f':
          fptr = (real*) (base + datum_dim_parms[i].off);
          *fptr = (real) atof(parm);
          return 1;

        /* string */
        case 's':
          sptr = (char*) (base + datum_dim_parms[i].off);
          memcpy(sptr, parm, 8);
          return 1;

        /* unknown */
        default:
          throw("invalid parameter type '%c'", datum_dim_parms[i].type);
      }
    }
  }

  /* no match was found. throw an exception. */
  throw("invalid dimension parameter name '%s'", name);
}

/* datum_dims_realloc(): reallocate the array of dimension parameters inside
 * a datum structure.
 * @D: pointer to the target datum structure.
 * @nd: new dimensionality
 */
int datum_dims_realloc (datum *D, unsigned int nd) {
  /* declare a loop counter. */
  unsigned int i;

  /* attempt to reallocate the dimensions array. */
  D->dims = (datum_dim*) realloc(D->dims, nd * sizeof(datum_dim));

  /* check that the reallocation succeeded. */
  if (D->dims == NULL)
    throw("failed to reallocate %u datum dimensions", nd);

  /* zero the contents of any new dimensions. */
  for (i = D->nd; i < nd; i++)
    memset(D->dims + i, 0, sizeof(datum_dim));

  /* set the new dimensionality. */
  D->nd = nd;

  /* return success. */
  return 1;
}

/* datum_dims_reorder(): reorder the dimensions inside a datum structure
 * according to the initial ordering @order. the target ordering will be
 * achieved by sorting the values in @order until they are in increasing
 * order.
 *
 * @D: pointer to the target datum structure.
 * @order: initial dimension ordering.
 */
int datum_dims_reorder (datum *D, int *order) {
  /* declare required variables:
   * @i: insertion sort outer loop counter.
   * @j: insertion sort inner loop counter.
   * @dim: a swap location for dimension information.
   */
  unsigned int i, j;
  datum_dim dim;
  int *ord, swp;

  /* allocate a temporary array of indices. */
  ord = (int*) malloc(D->nd * sizeof(int));
  if (!ord)
    throw("failed to allocate %u indices", D->nd);

  /* copy the ordering into the temporary array. */
  memcpy(ord, order, D->nd * sizeof(int));

  /* loop over the array of dimensions. */
  for (i = 1; i < D->nd; i++) {
    /* set the initial inner loop index. */
    j = i;

    /* loop over the unsorted dimensions. */
    while (j > 0 && ord[j - 1] > ord[j]) {
      /* swap the (j-1) and (j) elements. */
      dim = D->dims[j];
      D->dims[j] = D->dims[j - 1];
      D->dims[j - 1] = dim;

      /* swap the ordering of the indices. */
      swp = ord[j];
      ord[j] = ord[j - 1];
      ord[j - 1] = swp;

      /* decrement the inner loop counter. */
      j--;
    }
  }

  /* free the ordering array. */
  free(ord);

  /* return success. */
  return 1;
}

