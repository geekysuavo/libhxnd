
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

/* fn_def: function definition structure for holding all available function
 * names and pointers.
 */
struct fn_def {
  /* @name: the function name string.
   * @fn: the function pointer.
   */
  const char *name;
  int (*fn) (datum *D, const int dim, const char *argstr);
};

/* functions: local table of all available function names and pointers.
 */
static const struct fn_def functions[] = {
  { FN_NAME_ABS,      &fn_execute_abs },
  { FN_NAME_ADD,      &fn_execute_add },
  { FN_NAME_CUT,      &fn_execute_cut },
  { FN_NAME_FFT,      &fn_execute_fft },
  { FN_NAME_HT,       &fn_execute_ht },
  { FN_NAME_IST,      &fn_execute_ist },
  { FN_NAME_PHASE,    &fn_execute_phase },
  { FN_NAME_RESIZE,   &fn_execute_resize },
  { FN_NAME_SCALE,    &fn_execute_scale },
  { FN_NAME_SHIFT,    &fn_execute_shift },
  { FN_NAME_WINDOW,   &fn_execute_window },
  { FN_NAME_ZEROFILL, &fn_execute_zerofill },
  { NULL, NULL }
};

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

/* fn_parse_arg_intarray(): parses an array of integers from a key-value
 * string array.
 * @v: the string array.
 * @c: the string array length.
 */
int *fn_parse_arg_intarray (char **v, int c) {
  /* declare a few required variables:
   */
  unsigned int i, n;
  char **strv;
  int *intv;

  /* check the array length. */
  if (c != 2) {
    /* raise an exception and return nothing. */
    raise("r-value required for int-array parsing");
    return NULL;
  }

  /* convert the leading character to removable whitespace. */
  if (v[1][0] == '(' || v[1][0] == '{')
    v[1][0] = ' ';

  /* convert the trailing character to removable whitespace. */
  n = (strlen(v[1]) > 0 ? strlen(v[1]) - 1 : 0);
  if (v[1][n] == ')' || v[1][n] == '}')
    v[1][n] = ' ';

  /* split the value from the string array into a new string array. */
  strv = strsplit(v[1], ".", &n);

  /* check that the split was successful. */
  if (!strv || n < 1) {
    /* raise an exception and return nothing. */
    raise("failed to split int-array string");
    return NULL;
  }

  /* trim whitespace from the string array. */
  strvtrim(strv, n);

  /* allocate an array of integers having the same number of elements as
   * the string array.
   */
  intv = (int*) calloc(n + 1, sizeof(int));
  if (!intv) {
    /* raise an exception and return nothing. */
    raise("failed to allocate int-array");
    return NULL;
  }

  /* loop over the strings in the string array. */
  for (i = 0, intv[0] = n; i < n; i++)
    intv[i + 1] = atoi(strv[i]);

  /* free the string array */
  strvfree(strv, n);

  /* return the complete array of integers. */
  return intv;
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

  /* parse the value from the string array.
   */
  return atof(v[1]);
}

/* fn_parse_arg_str(): parses a string argument from a key-value string array.
 * this allocates a new string in memory, which needs to be freed by the
 * calling function of fn_scan_args().
 * @v: the string array.
 * @c: the string array length.
 */
char *fn_parse_arg_str (char **v, int c) {
  /* declare a few required variables. */
  unsigned int n;
  char *str;

  /* check the array length. */
  if (c != 2) {
    /* raise an error and return nothing. */
    raise("r-value required for string parsing");
    return NULL;
  }

  /* get the length of the input string. */
  n = strlen(v[1]) + 2;

  /* allocate memory for the string output. */
  str = (char*) malloc(n * sizeof(char));
  if (str)
    strcpy(str, v[1]);

  /* return the newly allocated, copied string. */
  return str;
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

      /* int array. */
      case FN_ARGTYPE_INTS:
        *((int**) ptr) = (found ? fn_parse_arg_intarray(valv, valc) : NULL);
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

      /* string. */
      case FN_ARGTYPE_STRING:
        *((char**) ptr) = (found ? fn_parse_arg_str(valv, valc) : NULL);
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
  /* declare a few required variables. */
  int i, len1, len2, len, matches, last_match;

  /* initialize the match counter and index. */
  matches = 0;
  last_match = 0;

  /* loop over the array of possible function names. */
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
    /* awesome. execute the function. */
    return functions[last_match].fn(D, dim, argstr);
  }
  else if (matches == 0) {
    /* no match. */
    throw("no functions matching name '%s'", name);
  }
  else if (matches > 1) {
    /* multiple matches. */
    throw("function name '%s' is ambiguous", name);
  }

  /* this should never occur. */
  throw("unknown error calling function '%s'", name);
}

