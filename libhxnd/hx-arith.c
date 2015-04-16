
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

/* * * * * * * * * * * RAW DATA OPERATIONS * * * * * * * * * * */

/* hx_data_add(): add the raw array elements of two hypercomplex values.
 * @xa: the raw array data of the first operand.
 * @xb: the raw array data of the second operand.
 * @xc: the raw array data of the result.
 * @s: the real factor to apply to b during addition.
 * @d: the algebraic dimensionality of the scalars.
 * @n: the number of array elements of the scalars.
 *
 * operation:
 *   c <= a + s * b
 *
 * operands:
 *   a: hypercomplex scalar.
 *   b: hypercomplex scalar.
 *   s: real scalar.
 *   c: hypercomplex scalar.
 */
int hx_data_add (real *xa, real *xb, real *xc, real s, int d, int n) {
  /* declare a required variable. */
  int i;

  /* check if a modified operation is to be performed. */
  if (xa == NULL) {
    /* modified: only scale b[i]. */
    for (i = 0; i < n; i++)
      xc[i] = s * xb[i];
  }
  else if (xb == NULL) {
    /* modified: only add s to a[i]. (i.e. b[i] == 1) */
    for (i = 0; i < n; i++)
      xc[i] = xa[i] + s;
  }
  else {
    /* full: scale the value of b[i] and add the result to a[i] */
    for (i = 0; i < n; i++)
      xc[i] = xa[i] + s * xb[i];
  }

  /* return success. */
  return 1;
}

/* hx_data_mul(): multiply the raw array elements of two hypercomplex values.
 * @xa: the raw array data of the first operand.
 * @xb: the raw array data of the second operand.
 * @xc: the raw array data of the result.
 * @d: the algebraic dimensionality of the scalars.
 * @n: the number of array elements of the scalars.
 * @tbl: the multiplication table to utilize.
 *
 * operation:
 *   c <= c + a * b
 *
 * operands:
 *   a: hypercomplex scalar.
 *   b: hypercomplex scalar.
 *   c: hypercomplex scalar.
 */
int hx_data_mul (real *xa, real *xb, real *xc, int d, int n, hx_algebra tbl) {
  /* declare a few required variables. */
  int i, j, k, tij;
  real sgn;

  /* loop over the first value's array elements. */
  for (i = 0; i < n; i++) {
    /* skip zero-valued left-operand coefficients. */
    if (xa[i] == 0.0)
      continue;

    /* loop over the second value's array elements. */
    for (j = 0; j < n; j++) {
      /* get the output index from the algebra table. */
      tij = tbl[i * n + j];
      k = (tij > 0 ? tij : -tij);

      /* get the output sign from the algebra table. */
      sgn = (real) (tij / k);
      k--;

      /* compute the final result. */
      xc[k] += sgn * xa[i] * xb[j];
    }
  }

  /* return success. */
  return 1;
}

/* hx_data_shuf(): shuffle the raw array elements of two hypercomplex values
 * in order to extract true hypercomplex information from existing gradient-
 * enhanced information.
 * @xa: the raw array data of the first input operand.
 * @xb: the raw array data of the second input operand.
 * @xc: the raw array data of the first output operand.
 * @xd: the raw array data of the second output operand.
 * @xph: the raw array data of a temporary phasor.
 *
 * operation:
 *   c <= (b + a)
 *   d <= (b - a) * u_k
 *
 * operands:
 *   a: hypercomplex scalar.
 *   b: hypercomplex scalar.
 *   c: hypercomplex scalar.
 *   d: hypercomplex scalar.
 *   ph: hypercomplex scalar.
 *
 * NOTE: this function may not be used to perform in-place shuffling.
 */
int hx_data_shuf (real *xa, real *xb, real *xc, real *xd,
                  real *xph, real *xtmp, int d, int n,
                  hx_algebra tbl) {
  /* compute the sum. */
  if (!hx_data_add(xb, xa, xc, 1.0, d, n))
    return 0;

  /* compute the difference. */
  if (!hx_data_add(xb, xa, xtmp, -1.0, d, n))
    return 0;

  /* explicitly compute the phasor array values. */
  hx_data_zero(xph, n);
  xph[1 << (d - 1)] = 1.0;

  /* phase-shift the difference. */
  hx_data_zero(xd, n);
  if (!hx_data_mul(xtmp, xph, xd, d, n, tbl))
    return 0;

  /* return success. */
  return 1;
}

/* hx_data_copy(): copy the coefficients between the raw arrays
 * of two hypercomplex values.
 * @x: the raw array data of the input operand.
 * @xcpy: the raw array data of the output operand.
 * @n: the number of array elements of the operands.
 */
int hx_data_copy (real *x, real *xcpy, int n) {
  /* copy the values. */
  memcpy(xcpy, x, n * sizeof(real));

  /* return success. */
  return 1;
}

/* hx_data_conj(): negate all the imaginary coefficients in the raw array
 * of a hypercomplex value.
 * @x: the raw array data of the input operand.
 * @xh: the raw array data of the output operand.
 * @n: the number of array elements of the operands.
 */
int hx_data_conj (real *x, real *xh, int n) {
  /* declare a required variable. */
  int i;

  /* loop over the imaginary elements. */
  for (i = 1; i < n; i++)
    xh[i] = -1.0 * x[i];

  /* return success. */
  return 1;
}

/* hx_data_semiconj(): negate all the imaginary basis elements in the raw
 * array of a hypercomplex value.
 * @x: the raw array data of the input operand.
 * @xh: the raw array data of the output operand.
 * @d: the algebraic dimensionality of the operands.
 * @n: the number of array elements of the operands.
 */
int hx_data_semiconj (real *x, real *xh, int d, int n) {
  /* declare a required variable. */
  int i, dneg, ineg;

  /* copy the data into the output array. */
  hx_data_copy(x, xh, n);

  /* loop over the basis elements. */
  for (dneg = 0; dneg < d; dneg++) {
    /* compute the negation index. */
    ineg = 1 << dneg;

    /* loop over the imaginary coefficients. */
    for (i = ineg; i < n; i++) {
      /* check if the current index needs negation. */
      if (i & ineg)
        xh[i] *= -1.0;
    }
  }

  /* return success. */
  return 1;
}

/* hx_data_zero(): sets the raw coefficients of a hypercomplex value to zero.
 * @x: the raw array data of the input operand.
 * @n: the number of array elements of the operand.
 */
int hx_data_zero (real *x, int n) {
  /* set the coefficients to zero. */
  memset(x, 0, n * sizeof(real));

  /* return success. */
  return 1;
}

/* hx_data_fill(): sets the raw coefficients of a hypercomplex value.
 * @x: the raw array data of the input operand.
 * @n: the number of array elements of the operand.
 * @val: the value to fill the array elements with.
 */
int hx_data_fill (real *x, int n, real val) {
  /* declare a required variable. */
  int i;

  /* loop over the array elements. */
  for (i = 0; i < n; i++)
    x[i] = val;

  /* return success. */
  return 1;
}

/* hx_data_norm(): compute the norm of the raw array elements of
 * a hypercomplex value.
 * @x: the raw array data of the input operand.
 * @n: the number of array elements of the operand.
 *
 * operation:
 *   x <= ||x||
 *
 * operands:
 *   x: hypercomplex scalar.
 */
int hx_data_norm (real *x, int n) {
  /* declare a few required variables. */
  int i;

  /* square the first array element. */
  x[0] *= x[0];

  /* loop over the remaining array elements. */
  for (i = 1; i < n; i++) {
    /* sum in the product and zero the element. */
    x[0] += x[i] * x[i];
    x[i] = 0.0;
  }

  /* square-root the first (real) array element. */
  x[0] = sqrt(x[0]);

  /* return success. */
  return 1;
}

/* hx_data_real_norm(): compute the norm of the raw array elements of a
 * hypercomplex value, and return the result as a real number.
 * @x: the raw array data of the input operand.
 * @n: the number of array elements of the operand.
 *
 * operation:
 *   nrm <= ||x||
 *
 * operands:
 *   x: hypercomplex scalar.
 */
real hx_data_real_norm (real *x, int n) {
  /* declare a few required variables. */
  real nrm;
  int i;

  /* loop over the array elements. */
  for (i = 0, nrm = 0.0; i < n; i++)
    nrm += x[i] * x[i];

  /* return the norm. */
  return sqrt(nrm);
}

/* hx_data_negate_basis(): negates the specified basis element of the raw
 * array elements of a hypercomplex value.
 * @x: the raw array data of the input operand.
 * @d: the algebraic dimensionality of the operand.
 * @n: the number of array elements of the operand.
 * @dneg: the algebraic dimension to negate.
 */
int hx_data_negate_basis (real *x, int d, int n, int dneg) {
  /* declare a required variable. */
  int i, ineg;

  /* compute the binary mask of the negation dimension. */
  ineg = 1 << dneg;

  /* loop over the array elements. */
  for (i = 0; i < n; i++) {
    /* check if the current index needs negation. */
    if (i & ineg)
      x[i] *= -1.0;
  }

  /* return success. */
  return 1;
}

/* hx_array_reorder_bases(): reorders the basis elements of each scalar value
 * in a hypercomplex array. the target ordering will be achieved by swapping
 * pairs of bases in @order until they are in increasing order.
 *
 * NOTE 1: the @d values in @order *MUST* be unique and lie in [0, @d).
 * NOTE 2: the values in @order get re-arranged by this function.
 *
 * @x: the structure pointer to the target array.
 * @order: the dimension ordering array.
 */
int hx_data_reorder_bases (real *x, int d, int n, int *order) {
  /* declare a few required variables:
   * @di: first basis element index in each swap.
   * @dj: second basis element index in each swap.
   */
  int di, dj, ni, nj, i, j, iswp;
  real swp;

  /* loop over the basis elements. */
  for (di = 0; di < d - 1; di++) {
    /* if the current basis element is already in it's proper place, proceed
     * to the next basis element.
     */
    if (order[di] == di)
      continue;

    /* locate the other basis element in the order array that will swap it's
     * position with the current basis element.
     */
    for (dj = di + 1; dj < d; dj++) {
      if (order[dj] == di)
        break;
    }

    /* compute the coefficient masks for each basis element in the swap. */
    ni = 1 << di;
    nj = 1 << dj;

    /* loop over the coefficients. */
    for (i = 0; i < n; i++) {
      /* check if we're at a swappable coefficient. */
      if ((i & ni) && !(i & nj)) {
        /* yes! compute the other swap coefficient index. */
        j = (i & ~ni) | nj;

        /* swap the coefficients. */
        swp = x[i];
        x[i] = x[j];
        x[j] = swp;
      }
    }

    /* swap the order values. */
    iswp = order[di];
    order[di] = order[dj];
    order[dj] = iswp;
  }

  /* return success. */
  return 1;
}

/* * * * * * * * * * * SCALAR OPERATIONS * * * * * * * * * * */

/* hx_scalar_add(): add two hypercomplex scalar values.
 * @a: the structure pointer of the first operand.
 * @b: the structure pointer of the second operand.
 * @s: the real factor to apply to b during addition.
 * @c: the structure pointer of the result.
 *
 * operation:
 *   c <= a + s * b
 *
 * operands:
 *   a: hypercomplex scalar.
 *   b: hypercomplex scalar.
 *   s: real scalar.
 *   c: hypercomplex scalar.
 */
int hx_scalar_add (hx_scalar *a, hx_scalar *b, real s, hx_scalar *c) {
  /* check if the algebraic dimensionalities match. */
  if (hx_scalar_dims_cmp(a, b) != 0 ||
      hx_scalar_dims_cmp(a, c) != 0)
    throw("scalar algebraic dimension mismatch");

  /* perform the raw data operation. */
  return hx_data_add(a->x, b->x, c->x, s, a->d, a->n);
}

/* hx_scalar_mul(): multiply two hypercomplex scalar values.
 * @a: the structure pointer of the first operand.
 * @b: the structure pointer of the second operand.
 * @c: the structure pointer of the result.
 *
 * operation:
 *   c <= c + a * b
 *
 * operands:
 *   a: hypercomplex scalar.
 *   b: hypercomplex scalar.
 *   c: hypercomplex scalar.
 */
int hx_scalar_mul (hx_scalar *a, hx_scalar *b, hx_scalar *c) {
  /* check if the algebraic dimensionalities match. */
  if (hx_scalar_dims_cmp(a, b) != 0 ||
      hx_scalar_dims_cmp(a, c) != 0)
    throw("scalar algebraic dimension mismatch");

  /* perform the raw data operation. */
  return hx_data_mul(a->x, b->x, c->x, a->d, a->n, a->tbl);
}

/* hx_scalar_scale(): multiply a hypercomplex scalar and a real value.
 * @a: the structure pointer of the first operand.
 * @s: the value of the second operand.
 * @b: the structure pointer of the result.
 *
 * operation:
 *   b <= a * s
 *
 * operands:
 *   a: hypercomplex scalar.
 *   s: real scalar.
 *   b: hypercomplex scalar.
 */
int hx_scalar_scale (hx_scalar *a, real s, hx_scalar *b) {
  /* check if the algebraic dimensionalities match. */
  if (hx_scalar_dims_cmp(a, b) != 0)
    throw("scalar algebraic dimension mismatch");

  /* perform the raw data operation. */
  return hx_data_add(NULL, a->x, b->x, s, a->d, a->n);
}

/* hx_scalar_zero(): sets the coefficients of a hypercomplex scalar to zero.
 * @a: the structure pointer of the input operand.
 */
int hx_scalar_zero (hx_scalar *a) {
  /* set the coefficients to zero. */
  memset(a->x, 0, a->n * sizeof(real));

  /* return success. */
  return 1;
}

/* hx_scalar_fill(): set all the coefficients of a hypercomplex scalar.
 * @a: the structure pointer of the input operand.
 * @val: the value to fill the operand with.
 */
int hx_scalar_fill (hx_scalar *a, real val) {
  /* perform the raw data operation. */
  return hx_data_fill(a->x, a->n, val);
}

/* hx_scalar_norm(): compute the norm of a hypercomplex scalar.
 * @a: the structure pointer of the input operand.
 *
 * operation:
 *   a <= ||a||
 *
 * operands:
 *   a: hypercomplex scalar.
 */
int hx_scalar_norm (hx_scalar *a) {
  /* perform the raw data operation. */
  return hx_data_norm(a->x, a->n);
}

/* hx_scalar_negate_basis(): negates the specified basis element of a
 * hypercomplex scalar value.
 * @x: the structure pointer of the input operand.
 * @dneg: the algebraic dimension to negate.
 */
int hx_scalar_negate_basis (hx_scalar *x, int dneg) {
  /* check that the algebraic dimension is in bounds. */
  if (dneg < 0 || dneg >= x->d)
    throw("algebraic dimension %d out of bounds [0,%d)", dneg, x->d);

  /* perform the raw data operation. */
  return hx_data_negate_basis(x->x, x->d, x->n, dneg);
}

/* hx_scalar_reorder_bases(): reorders the basis elements of a hypercomplex
 * scalar value.
 * @x: the structure pointer to the target scalar.
 * @order: the dimension ordering array.
 */
int hx_scalar_reorder_bases (hx_scalar *x, int *order) {
  /* declare a few required variables. */
  int *scratch;
  int i, ret;

  /* allocate scratch space for the reordering operation. */
  scratch = (int*) malloc(x->d * sizeof(int));
  if (!scratch)
    throw("failed to allocate scratch space");

  /* check the bounds on the ordering array. */
  for (i = 0; i < x->d; i++) {
    /* check that the current order is in bounds. */
    if (order[i] < 0 || order[i] >= x->d)
      throw("order %d (#%d) out of bounds [0,%d)", order[i], i, x->d);
  }

  /* copy the current ordering into the scratch space. */
  memcpy(scratch, order, x->d * sizeof(int));

  /* perform the raw scalar data operation. */
  ret = hx_data_reorder_bases(x->x, x->d, x->n, scratch);

  /* free the scratch space and return the result. */
  free(scratch);
  return ret;
}

/* * * * * * * * * * * ARRAY OPERATIONS * * * * * * * * * * */

/* hx_array_add_scalar(): add a hypercomplex array and a scalar.
 * @a: the structure pointer to the array operand.
 * @b: the structure pointer to the scalar operand.
 * @s: the real factor to apply to b during addition.
 * @c: the structure pointer to the result.
 *
 * operation:
 *   c <= a + s * b
 *
 * operands:
 *   a: hypercomplex array.
 *   b: hypercomplex scalar.
 *   s: real scalar.
 *   c: hypercomplex array.
 */
int hx_array_add_scalar (hx_array *a, hx_scalar *b, real s, hx_array *c) {
  /* declare a required variable. */
  int i;

  /* check if the algebraic dimensionalities match. */
  if (a->d != b->d || hx_array_conf_cmp(a, c) != 0)
    throw("array-scalar configuration mismatch");

  /* loop over the array elements. */
  for (i = 0; i < a->len; i += a->n) {
    /* perform the raw scalar data operation. */
    if (!hx_data_add(a->x + i, b->x, c->x + i, s, a->d, a->n))
      return 0;
  }

  /* return success. */
  return 1;
}

/* hx_array_add_array(): add two hypercomplex arrays.
 * @a: the structure pointer to the first operand.
 * @b: the structure pointer to the second operand.
 * @s: the real factor to apply to b during addition.
 * @c: the structure pointer to the result.
 *
 * operation:
 *   c <= a + s * b
 *
 * operands:
 *   a: hypercomplex array.
 *   b: hypercomplex array.
 *   s: real scalar.
 *   c: hypercomplex array.
 */
int hx_array_add_array (hx_array *a, hx_array *b, real s, hx_array *c) {
  /* declare a required variable. */
  int i;

  /* check if the array configurations match. */
  if (hx_array_conf_cmp(a, b) != 0 ||
      hx_array_conf_cmp(a, c) != 0)
    throw("array configuration mismatch");

  /* loop over the array elements. */
  for (i = 0; i < a->len; i += a->n) {
    /* perform the raw scalar data operation. */
    if (!hx_data_add(a->x + i, b->x + i, c->x + i, s, a->d, a->n))
      return 0;
  }

  /* return success. */
  return 1;
}

/* hx_array_mul_scalar(): multiply a hypercomplex array and a scalar.
 * @a: the structure pointer to the array operand.
 * @b: the structure pointer to the scalar operand.
 * @c: the structure pointer to the result.
 *
 * operation:
 *   c <= c + a * b
 *
 * operands:
 *   a: hypercomplex array.
 *   b: hypercomplex scalar.
 *   c: hypercomplex array.
 */
int hx_array_mul_scalar (hx_array *a, hx_scalar *b, hx_array *c) {
  /* declare a required variable. */
  int i;

  /* check if the algebraic dimensionalities match. */
  if (a->d != b->d || hx_array_conf_cmp(a, c) != 0)
    throw("array-scalar configuration mismatch");

  /* loop over the array elements. */
  for (i = 0; i < a->len; i += a->n) {
    /* perform the raw scalar data operation. */
    if (!hx_data_mul(b->x, a->x + i, c->x + i, a->d, a->n, a->tbl))
      return 0;
  }

  /* return success. */
  return 1;
}

/* hx_array_mul_array(): multiply two hypercomplex arrays.
 * @a: the structure pointer to the first operand.
 * @b: the structure pointer to the second operand.
 * @c: the structure pointer to the result.
 *
 * operation:
 *   c <= c + a * b
 *
 * operands:
 *   a: hypercomplex array.
 *   b: hypercomplex array.
 *   c: hypercomplex array.
 */
int hx_array_mul_array (hx_array *a, hx_array *b, hx_array *c) {
  /* declare a required variable. */
  int i;

  /* check if the array configurations match. */
  if (hx_array_conf_cmp(a, b) != 0 ||
      hx_array_conf_cmp(a, c) != 0)
    throw("array configuration mismatch");

  /* loop over the array elements. */
  for (i = 0; i < a->len; i += a->n) {
    /* perform the raw scalar data operation. */
    if (!hx_data_mul(a->x + i, b->x + i, c->x + i, a->d, a->n, a->tbl))
      return 0;
  }

  /* return success. */
  return 1;
}

/* hx_array_mul_vector_cb(): callback function for hx_array_mul_vector().
 *
 * args:
 *  see hx_array_foreach_cb().
 *
 * varargs:
 *  @b: vector multiplier.
 *  @ytmp: temporary array, same configuration as @b.
 */
int hx_array_mul_vector_cb (hx_array *x, hx_array *y,
                            hx_index idx, int pidx,
                            va_list *vl) {
  /* extract the varargs. */
  hx_array *b = va_arg(*vl, hx_array*);
  hx_array *ytmp = va_arg(*vl, hx_array*);

  /* copy the array contents from @y to @ytmp. */
  memcpy(ytmp->x, y->x, y->len * sizeof(real));

  /* zero the destination array coefficients. */
  hx_array_zero(y);

  /* perform the array multiplication. */
  if (!hx_array_mul_array(ytmp, b, y))
    throw("failed to multiply vector inside array");

  /* return success. */
  return 1;
}

/* hx_array_mul_vector(): multiply each vector of an array @a along a given
 * (topological) dimension @kmul by a vector-shaped array @b.
 *
 * @a: the structure pointer to the first (whole-array) operand.
 * @b: the structure pointer to the second (vector) operand.
 * @kmul: the dimension index along which to multiply.
 * @c: the structure pointer to the result.
 *
 * operation:
 *   c_k <= a_k * b, foreach k in a.sz[kmul]
 *
 * operands:
 *   a: hypercomplex array, at least one-dimensional.
 *   b: hypercomplex array, must be one-dimensional.
 *   c: hypercomplex array, same dimensionality as @a.
 *
 * NOTE: the array pointers @a and @c are permitted to be the same value,
 *       effectively allowing in-place multiplication of @a by @b.
 */
int hx_array_mul_vector (hx_array *a, hx_array *b, int kmul, hx_array *c) {
  /* declare a required variable:
   * @ytmp: temporary array of same configuration as slices.
   */
  hx_array ytmp;

  /* check that the dimension index is in bounds. */
  if (kmul < 0 || kmul >= a->k)
    throw("array index %d out of bounds [0,%d)", kmul, a->k);

  /* check that the vector operand is truly a vector. */
  if (!hx_array_is_vector(b))
    throw("vector argument has invalid dimensionality");

  /* check that the array algebraic dimensionalities match. */
  if (hx_array_dims_cmp(a, b) ||
      hx_array_dims_cmp(a, c))
    throw("array algebraic dimensionality mismatch");

  /* check that the important topological dimensionalities match. */
  if (a->sz[kmul] != b->sz[0] ||
      hx_array_topo_cmp(a, c))
    throw("array topological dimensionality mismatch");

  /* allocate a temporary array for destination slice storage. */
  if (!hx_array_copy(&ytmp, b))
    throw("failed to allocate slice (%d, 1)-array", b->d);

  /* copy the array contents from the first operand to the destination. */
  if (a != c)
    memcpy(c->x, a->x, a->len * sizeof(real));

  /* run the callback function over every vector along @kmul. */
  if (!hx_array_foreach_vector(c, kmul, &hx_array_mul_vector_cb, b, &ytmp))
    return 0;

  /* free the temporary array. */
  hx_array_free(&ytmp);

  /* return success. */
  return 1;
}

/* hx_array_scale(): multiply a hypercomplex array and a real value.
 * @a: the structure pointer of the first operand.
 * @s: the value of the second operand.
 * @b: the structure pointer of the result.
 *
 * operation:
 *   b <= a * s
 *
 * operands:
 *   a: hypercomplex array.
 *   s: real scalar.
 *   b: hypercomplex array.
 */
int hx_array_scale (hx_array *a, real s, hx_array *b) {
  /* declare a required variable. */
  int i;

  /* check if the array configurations match. */
  if (hx_array_conf_cmp(a, b) != 0)
    throw("array configuration mismatch");

  /* loop over the array elements. */
  for (i = 0; i < a->len; i += a->n) {
    /* perform the raw scalar data operation. */
    if (!hx_data_add(NULL, a->x + i, b->x + i, s, a->d, a->n))
      return 0;
  }

  /* return success. */
  return 1;
}

/* hx_array_zero(): sets the coefficients of a hypercomplex array to zero.
 * @a: the structure pointer of the input operand.
 */
int hx_array_zero (hx_array *a) {
  /* set the coefficients to zero. */
  memset(a->x, 0, a->len * sizeof(real));

  /* return success. */
  return 1;
}

/* hx_array_fill(): set all the coefficients of a hypercomplex array.
 * @a: the structure pointer of the input operand.
 * @val: the value to fill the operand with.
 */
int hx_array_fill (hx_array *a, real val) {
  /* perform the raw data operation. */
  return hx_data_fill(a->x, a->len, val);
}

/* hx_array_norm(): compute the norm every scalar in a hypercomplex array.
 * @a: the structure pointer to the input operand.
 *
 * operation:
 *   a <= ||a||
 *
 * operands:
 *   a: hypercomplex array.
 */
int hx_array_norm (hx_array *a) {
  /* declare a required variable. */
  int i;

  /* loop over the array elements. */
  for (i = 0; i < a->len; i += a->n) {
    /* perform the raw scalar data operation. */
    if (!hx_data_norm(a->x + i, a->n))
      return 0;
  }

  /* return success. */
  return 1;
}

/* hx_array_alternate_sign_cb(): callback function for
 * hx_array_alternate_sign().
 *
 * args:
 *  see hx_array_foreach_cb().
 *
 * varargs:
 *  none.
 */
int hx_array_alternate_sign_cb (hx_array *x, hx_array *y,
                                int *arr, int idx,
                                va_list *vl) {
  /* declare a required variable. */
  int i;

  /* loop over the array elements. */
  for (i = y->n; i < y->len; i += 2 * y->n) {
    /* perform the raw scalar data operation. */
    if (!hx_data_add(NULL, y->x + i, y->x + i, -1.0, y->d, y->n))
      return 0;
  }

  /* return success. */
  return 1;
}

/* hx_array_alternate_sign(): negate every other scalar in a hypercomplex
 * array, in place, along a specified topological dimension.
 * @x: the structure pointer to the input operand.
 * @k: dimension along which alternation will be performed.
 */
int hx_array_alternate_sign (hx_array *x, int k) {
  /* check that the shift dimension is in bounds. */
  if (k < 0 || k >= x->k)
    throw("alternation dimension %d out of bounds [0,%d)", k, x->k);

  /* perform the per-vector sign alternation operation. */
  if (!hx_array_foreach_vector(x, k, &hx_array_alternate_sign_cb))
    throw("failed to perform alternation");

  /* return success. */
  return 1;
}

/* hx_array_negate_basis(): negates the specified basis element of each
 * scalar value in a hypercomplex array.
 * @x: the structure pointer to the target array.
 * @dneg: the algebraic dimension to negate.
 */
int hx_array_negate_basis (hx_array *x, int dneg) {
  /* declare a required variable. */
  int i;

  /* check that the algebraic dimension is in bounds. */
  if (dneg < 0 || dneg >= x->d)
    throw("algebraic dimension %d out of bounds [0,%d)", dneg, x->d);

  /* loop over the array elements. */
  for (i = 0; i < x->len; i += x->n) {
    /* perform the raw scalar data operation. */
    if (!hx_data_negate_basis(x->x + i, x->d, x->n, dneg))
      return 0;
  }

  /* return success. */
  return 1;
}

/* hx_array_reorder_bases(): reorders the basis elements of each scalar value
 * in a hypercomplex array.
 * @x: the structure pointer to the target array.
 * @order: the dimension ordering array.
 */
int hx_array_reorder_bases (hx_array *x, int *order) {
  /* declare a required variable. */
  int *scratch;
  int i;

  /* allocate scratch space for the reordering operation. */
  scratch = (int*) malloc(x->d * sizeof(int));
  if (!scratch)
    throw("failed to allocate scratch space");

  /* check the bounds on the ordering array. */
  for (i = 0; i < x->d; i++) {
    /* check that the current order is in bounds. */
    if (order[i] < 0 || order[i] >= x->d)
      throw("order %d (#%d) out of bounds [0,%d)", order[i], i, x->d);
  }

  /* loop over the array elements. */
  for (i = 0; i < x->len; i += x->n) {
    /* copy the current ordering into the scratch space. */
    memcpy(scratch, order, x->d * sizeof(int));

    /* perform the raw scalar data operation. */
    if (!hx_data_reorder_bases(x->x + i, x->d, x->n, scratch))
      return 0;
  }

  /* free the scratch space and return success. */
  free(scratch);
  return 1;
}

