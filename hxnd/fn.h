
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

/* fn_valtype: enumerated type for all accepted function argument value types.
 */
enum fn_valtype {
  FN_VALTYPE_UNKNOWN = 0,
  FN_VALTYPE_INT     = 1,  /* integer:     @i  */
  FN_VALTYPE_INTS    = 2,  /* int-array:   @iv */
  FN_VALTYPE_BOOL    = 3,  /* boolean:     @b  */
  FN_VALTYPE_FLOAT   = 4,  /* float:       @f  */
  FN_VALTYPE_FLOATS  = 5,  /* float-array: @fv */
  FN_VALTYPE_STRING  = 6,  /* string:      @s  */
  FN_VALTYPE_CHUNK   = 7   /* chunk:       @p  */
};

/* fn_val: union of all values that function arguments may hold.
 */
union fn_val {
  int i, *iv, b;
  real f, *fv;
  char *s;
  void *p;
};

/* fn_arg: structure holding a single function argument, useful for defining
 * and lexing function arguments to and from strings.
 */
typedef struct {
  /* @name: argument name string.
   * @val: argument value union.
   * @sz: value size, in bytes.
   */
  const char *name;
  union fn_val val;
  size_t sz;

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

/* function declarations (fn.c): */

int fn_execute (void *fndata, const int dim, fn *func, fn_arg *args);

int fn_execute_from_strings (void *fndata, const int dim,
                             const char *fnname,
                             const char *argstr);

/* function declarations (fn-args.c): */

int fn_args_get (const fn_arg *argdef, const int i, void *val, size_t *sz);

int fn_args_set (fn_arg *argdef, const int i, void *val, size_t sz);

int fn_args_get_all (const fn_arg *argdef, ...);

int fn_args_get_sizes (const fn_arg *argdef, ...);

int fn_args_from_string (fn_arg *argdef, const char *argstr);

fn_arg *fn_args_copy (const fn_arg *argsrc);

/* function declarations (fn-list.c): */

void fn_list_init (fn_list *fl);

void fn_list_free (fn_list *fl);

#endif /* __HXND_FN_H__ */

