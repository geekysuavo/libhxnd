
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

/* fn_parse_arg_int(): parses an integer from a key-value string array.
 * @v: the string array.
 * @c: the string array length.
 */
int fn_parse_arg_int (char **v, int c) {
  /* check the array length. */
  if (c != 2)
    throw("r-value required for integer parsing");

  /* parse the value from the string array. */
  return atoi(v[1]);
}

/* fn_parse_arg_bool(): parses a boolean (in the form of an integer) from
 * a key-value string array.
 * @v: the string array.
 * @c: the string array length.
 */
int fn_parse_arg_bool (char **v, int c) {
  /* check the array length. */
  if (c == 1) {
    /* placing a boolean without an r-value in the arguments shall flag a
     * 'true' value, unless prefixed with 'no'.
     */
    if (strlen(v[0]) >= 3 && strncmp(v[0], "no", 2) == 0)
      return 0;
    else
      return 1;
  }
  else if (c == 2) {
    /* check for acceptable values for 'true'. */
    if (strcmp(v[1], "true") == 0 ||
        strcmp(v[1], "True") == 0 ||
        strcmp(v[1], "TRUE") == 0 ||
        strcmp(v[1], "yes") == 0 ||
        strcmp(v[1], "Yes") == 0 ||
        strcmp(v[1], "YES") == 0 ||
        strcmp(v[1], "y") == 0 ||
        strcmp(v[1], "Y") == 0)
      return 1;
    else
      return 0;
  }
  else
    throw("unsupported argument syntax for boolean parsing");

  /* this will never occur. */
  return 0;
}

/* fn_parse_arg_float(): parses a floating-point number from a key-value
 * string array.
 * @v: the string array.
 * @c: the string array length.
 */
real fn_parse_arg_float (char **v, int c) {
  /* check the array length. */
  if (c != 2)
    throw("r-value required for float parsing");

  /* parse the value from the string array. */
  return atof(v[1]);
}

/* fn_scan_args(): scans a processing function string and extracts all known
 * function arguments. throws an error when unknown arguments or invalid
 * argument values are supplied.
 * @argstr: the processing function string to scan.
 * @argdef: the argument definition structure. must be null-terminated.
 * @...: pointers to accepted arguments, in the order their definitions
 *       appear in @argdef.
 */
int fn_scan_args (const char *argstr, const fn_args *argdef, ...) {
  /* declare variables for variable argument handling:
   * @argi: the argument index.
   * @vl: the variable arguments list.
   */
  unsigned int argi;
  va_list vl;

  /* declare variables for argument list parsing:
   * @argv: the argument string array.
   * @argc: the argument string array length.
   */
  unsigned int argc;
  char **argv;

  /* declare variables for individual argument parsing:
   * @valv: argument key-value string array.
   * @valc: argument key-value string array length.
   * @defi: argument definition structure loop index.
   */
  unsigned int valc, defi, found;
  char **valv, *noname;
  void *ptr;

  /* parse the argument string into an array. */
  argv = strsplit(argstr, ",", &argc);

  /* check that the split was successful. */
  if (!argv || !argc)
    throw("failed to split argument string");

  /* trim the string array elements. */
  strvtrim(argv, argc);

  /* initialize the variable arguments list. */
  va_start(vl, argdef);
  valv = NULL;

  /* loop over the argument definition structure entries. */
  for (defi = 0; argdef[defi].name; defi++) {
    /* initialize the negated argument name. */
    noname = NULL;

    /* check if a negated argument should be searched for as well. */
    if (argdef[defi].type == FN_ARGTYPE_BOOL) {
      /* build the negated string. */
      noname = (char*) malloc((strlen(argdef[defi].name) + 4) * sizeof(char));
      if (noname)
        sprintf(noname, "no%s", argdef[defi].name);
    }

    /* search the elements of the argument string array. */
    for (argi = 0, found = 0; argi < argc; argi++) {
      /* parse the string into an l-value and an r-value. */
      valv = strsplit(argv[argi], "=", &valc);

      /* check that the split was successful. */
      if (!valv || valc < 1)
        throw("failed to split argument '%s'", argv[argi]);

      /* trim the string array elements. */
      strvtrim(valv, valc);

      /* compare the l-value against the current argument name. */
      if (strcmp(valv[0], argdef[defi].name) == 0 ||
          (noname && strlen(valv[0]) >= 3 &&
           strcmp(valv[0], noname) == 0)) {
        /* break the search. */
        found = 1;
        break;
      }

      /* free the values string array. */
      strvfree(valv, valc);
    }

    /* the search was successful. pull another pointer out of the
     * variable arguments list for result storage.
     */
    ptr = va_arg(vl, void*);

    /* act based on the expected type of the argument. */
    switch (argdef[defi].type) {
      /* signed int. */
      case FN_ARGTYPE_INT:
        *((int*) ptr) = (found ? fn_parse_arg_int(valv, valc) :
                                 atoi(argdef[defi].def));
        break;

      /* boolean. */
      case FN_ARGTYPE_BOOL:
        *((int*) ptr) = (found ? fn_parse_arg_bool(valv, valc) :
                                 atoi(argdef[defi].def));
        break;

      /* float. */
      case FN_ARGTYPE_FLOAT:
        *((real*) ptr) = (found ? fn_parse_arg_float(valv, valc) :
                                  atof(argdef[defi].def));
        break;

      /* other... */
      default:
        throw("unsupported argument type '%c'", argdef[defi].type);
    }

    /* free the values string array, if required. */
    if (found)
      strvfree(valv, valc);

    /* free the negated value string, if required. */
    if (noname)
      free(noname);
  }

  /* free the variable arguments list and the argument string array. */
  strvfree(argv, argc);
  va_end(vl);

  /* return success. */
  return 1;
}

/* fn_execute(): global processing function executor. runs the appropriate
 * fn_execute_*() processing function on the datum pointer @D, which will
 * be modified in-place.
 */
int fn_execute (datum *D,
                const char *name,
                const int dim,
                const char *argstr) {
  /* act based on the function name. */
  if (strcmp(name, FN_NAME_FFT) == 0) {
    /* execute the 'fft' function. */
    return fn_execute_fft(D, dim, argstr);
  }
  else
    throw("invalid function name '%s'", name);

  /* return success. */
  return 1;
}

