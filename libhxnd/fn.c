
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

/* include the processing function header. */
#include <hxnd/fn.h>

/* include the function argument definitions header. */
#include <hxnd/fn-args.h>

/* include the function pointer declarations header. */
#include <hxnd/fn-handlers.h>

/* functions: local table of all available function names, pointers
 * and argument definitions.
 */
static fn functions[] = {
  { FN_NAME_ABS,      (fn_pointer) &fn_abs,      NULL },
  { FN_NAME_ADD,      (fn_pointer) &fn_add,      fn_args_add },
  { FN_NAME_COMPLEX,  (fn_pointer) &fn_complex,  NULL },
  { FN_NAME_CUT,      (fn_pointer) &fn_cut,      fn_args_cut },
  { FN_NAME_FFT,      (fn_pointer) &fn_fft,      fn_args_fft },
  { FN_NAME_HT,       (fn_pointer) &fn_ht,       NULL },
  { FN_NAME_IST,      (fn_pointer) &fn_ist,      fn_args_ist },
  { FN_NAME_MULTIPLY, (fn_pointer) &fn_multiply, fn_args_multiply },
  { FN_NAME_PHASE,    (fn_pointer) &fn_phase,    fn_args_phase },
  { FN_NAME_REAL,     (fn_pointer) &fn_real,     NULL },
  { FN_NAME_RESIZE,   (fn_pointer) &fn_resize,   fn_args_resize },
  { FN_NAME_SHIFT,    (fn_pointer) &fn_shift,    fn_args_shift },
  { FN_NAME_WINDOW,   (fn_pointer) &fn_window,   fn_args_window },
  { FN_NAME_ZEROFILL, (fn_pointer) &fn_zerofill, fn_args_zerofill },
  { NULL,             NULL,         NULL }
};

/* fn_lookup(): execute a name-lookup of a processing function based on
 * a complete or partial string match.
 * @name: function name string to search with.
 */
fn *fn_lookup (const char *name) {
  /* declare a few required variables. */
  int i, len1, len2, len, matches, last_match;
  fn *func;

  /* initialize the match counter and index. */
  matches = 0;
  last_match = 0;

  /* initialize the function address. */
  func = NULL;

  /* loop over the array of available function names. */
  for (i = 0; functions[i].name; i++) {
    /* get the length of the strings in the comparison. */
    len1 = strlen(name);
    len2 = strlen(functions[i].name);

    /* get the minimum of the two lengths. */
    len = (len1 < len2 ? len1 : len2);

    /* check if the current function name is a match. */
    if (strncmp(name, functions[i].name, len) == 0) {
      /* yes. increment the match counter and store the index. */
      matches++;
      last_match = i;
    }
  }

  /* check if a match was identified. */
  if (matches == 1) {
    /* set the function address to the last matched address. */
    func = &functions[last_match];
  }
  else if (matches == 0) {
    /* no match. */
    raise("no functions matching name '%s'", name);
  }
  else if (matches > 1) {
    /* multiple matches. */
    raise("function name '%s' is ambiguous", name);
  }

  /* return the identified function. */
  return func;
}

/* fn_execute(): execute a processing function using a function pointer
 * and an argument definition array.
 * @fndata: pointer to the datum/dataset structure to manipulate in-place.
 * @dim: dimension to apply the function along, or -1.
 * @func: function structure pointer.
 * @args: function arguments.
 */
int fn_execute (void *fndata, const int dim, fn *func, fn_arg *args) {
  /* execute the function. pretty easy, eh? */
  return func->ptr(fndata, dim, args);
}

/* fn_execute_from_strings(): execute a processing function using string
 * representations of both the function name and arguments.
 * @fndata: pointer to the datum/dataset structure to manipulate in-place.
 * @dim: dimension to apply the function along, or -1.
 * @argstr: string of function arguments to parse.
 */
int fn_execute_from_strings (void *fndata, const int dim,
                             const char *fnname,
                             const char *argstr) {
  /* declare variables to hold parsed function arguments. */
  fn_arg *args = NULL;
  fn *func = NULL;
  int ret;

  /* look up the function address by its name. */
  func = fn_lookup(fnname);

  /* check if a match was identified. */
  if (!func)
    throw("failed to look up function '%s'", fnname);

  /* allocate a copy of the function's argdef array. */
  args = fn_args_copy(func->args);

  /* parse the argument string into the temporary argdef array. */
  if (!fn_args_from_string(args, argstr)) {
    /* free the allocated argdef array. */
    free(args);

    /* throw an exception. */
    throw("failed to parse argument string for '%s'", fnname);
  }

  /* execute the function. */
  ret = fn_execute(fndata, dim, func, args);

  /* free the allocated argdef array. */
  if (args)
    free(args);

  /* return the result of the function call. */
  return ret;
}

