
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

/* ensure once-only inclusion. */
#ifndef __HXND_FN_H__
#define __HXND_FN_H__

/* include the hypercomplex math headers. */
#include <hxnd/hx.h>

/* define string names for all available processing functions.
 */
#define FN_NAME_ABS       "abs"
#define FN_NAME_ADD       "add"
#define FN_NAME_COMPLEX   "complex"
#define FN_NAME_CUT       "cut"
#define FN_NAME_FFT       "fft"
#define FN_NAME_HT        "ht"
#define FN_NAME_IST       "ist"
#define FN_NAME_MULTIPLY  "multiply"
#define FN_NAME_PHASE     "phase"
#define FN_NAME_REAL      "real"
#define FN_NAME_RESIZE    "resize"
#define FN_NAME_SHIFT     "shift"
#define FN_NAME_WINDOW    "window"
#define FN_NAME_ZEROFILL  "zerofill"

/* define string names for all available multivariate functions.
 */
#define FN_NAME_NORM   "norm"
#define FN_NAME_SCALE  "scale"

/* fn_valtype: enumerated type for all accepted function argument value types.
 */
enum fn_valtype {
  FN_VALTYPE_UNKNOWN = 0,
  FN_VALTYPE_INT     = 1,  /* integer:     @i  */
  FN_VALTYPE_INTS    = 2,  /* int-array:   @iv */
  FN_VALTYPE_BOOL    = 3,  /* boolean:     @b  */
  FN_VALTYPE_FLOAT   = 4,  /* float:       @f  */
  FN_VALTYPE_FLOATS  = 5,  /* float-array: @fv */
  FN_VALTYPE_STRING  = 6   /* string:      @s  */
};

/* fn_val: union of all values that function arguments may hold.
 */
union fn_val {
  int i, *iv, b;
  real f, *fv;
  char *s;
};

/* fn_arg: structure holding a single function argument, useful for defining
 * and lexing function arguments to and from strings.
 */
typedef struct {
  /* @name: argument name string.
   * @val: argument value union.
   */
  const char *name;
  union fn_val val;

  /* @type: argument type. */
  enum fn_valtype type;
}
fn_arg;

/* fn_pointer: callback function prototype for all macro-executable
 * processing functions.
 * @fndata: pointer to the datum or dataset structure to manipulate.
 * @dim: dimensional to apply the function along, or -1.
 * @args: array of argument definitions and values.
 */
typedef int (*fn_pointer) (void *fndata, const int dim,
                           const fn_arg *args);

/* fn: structure holding a function pointer and an array of argument
 * definitions, i.e. a complete processing function definition.
 */
typedef struct {
  /* @name: function name string used for function lookup and display.
   * @ptr: function pointer to call when applying the function.
   * @args: array of argument definition structures.
   */
  const char *name;
  fn_pointer ptr;
  fn_arg *args;
}
fn;

/* fn_list: structure holding an array of functions, used as a means of
 * describing the total data handling path from raw data to processed
 * spectra, or from processed spectra to data matrices, etc.
 */
typedef struct {
  /* @v: array of function structures.
   * @n: number of elements in the array.
   */
  fn *v;
  int n;
}
fn_list;

/* include the nmr datum header. */
#include <hxnd/nmr-datum.h>

/* include the nmr dataset header. */
#include <hxnd/mx-dataset.h>

/* function declarations (fn.c): */

int fn_execute (void *fndata, const int dim, fn *func, fn_arg *args);

int fn_execute_from_strings (void *fndata, const int dim,
                             const char *fnname,
                             const char *argstr);

/* function declarations (fn-args.c): */

int fn_args_get (const fn_arg *argdef, const int i, void *val);

int fn_args_set (fn_arg *argdef, const int i, void *val);

int fn_args_get_all (const fn_arg *argdef, ...);

int fn_args_from_string (fn_arg *argdef, const char *argstr);

fn_arg *fn_args_copy (const fn_arg *argsrc);

/* function declarations (fn-list.c): */

void fn_list_init (fn_list *fl);

void fn_list_free (fn_list *fl);

/* function declarations, datum handlers: */

int fn_abs (datum *D, const int dim, const fn_arg *args);

int fn_add (datum *D, const int dim, const fn_arg *args);

int fn_complex (datum *D, const int dim, const fn_arg *args);

int fn_cut (datum *D, const int dim, const fn_arg *args);

int fn_fft (datum *D, const int dim, const fn_arg *args);

int fn_ht (datum *D, const int dim, const fn_arg *args);

int fn_ist (datum *D, const int dim, const fn_arg *args);

int fn_multiply (datum *D, const int dim, const fn_arg *args);

int fn_phase (datum *D, const int dim, const fn_arg *args);

int fn_real (datum *D, const int dim, const fn_arg *args);

int fn_resize (datum *D, const int dim, const fn_arg *args);

int fn_shift (datum *D, const int dim, const fn_arg *args);

int fn_window (datum *D, const int dim, const fn_arg *args);

int fn_zerofill (datum *D, const int dim, const fn_arg *args);

/* function declarations, dataset handlers: */

int fn_norm (dataset *Dset, const int dim, const fn_arg *args);

int fn_scale (dataset *Dset, const int dim, const fn_arg *args);

#endif /* __HXND_FN_H__ */

