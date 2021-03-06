
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

/* include the string functions header. */
#include <hxnd/str.h>

/* fn_args_get(): retrieve an argument value from an argument definition
 * array.
 * @argdef: the argument definition array.
 * @i: the argument index to retrieve.
 * @val: pointer to the output value.
 * @sz: pointer to the output size.
 */
int fn_args_get (const fn_arg *argdef, const int i, void *val, size_t *sz) {
  /* declare a required variable. */
  int n;

  /* ensure the argdef array index is in bounds. */
  for (n = 0; argdef[n].name; n++) {}
  if (i >= n)
    throw("argument index %d out of bounds [0,%d)", i, n);

  /* determine the type of the argument value. */
  switch (argdef[i].type) {
    /* integer. */
    case FN_VALTYPE_INT:
      *((int*) val) = argdef[i].val.i;
      break;

    /* int-array. */
    case FN_VALTYPE_INTS:
      *((int**) val) = argdef[i].val.iv;
      break;

    /* boolean. */
    case FN_VALTYPE_BOOL:
      *((int*) val) = argdef[i].val.b;
      break;

    /* float. */
    case FN_VALTYPE_FLOAT:
      *((real*) val) = argdef[i].val.f;
      break;

    /* float-array. */
    case FN_VALTYPE_FLOATS:
      *((real**) val) = argdef[i].val.fv;
      break;

    /* string. */
    case FN_VALTYPE_STRING:
      *((char**) val) = argdef[i].val.s;
      break;

    /* chunk. */
    case FN_VALTYPE_CHUNK:
      *((void**) val) = argdef[i].val.p;
      break;

    /* other. */
    default:
      throw("unsupported argument type");
  }

  /* return the argument size, if requested. */
  if (sz)
    *sz = argdef[i].sz;

  /* return success. */
  return 1;
}

/* fn_args_set(): store an argument value into an argument definition array.
 * @argdef: the argument definition array.
 * @i: the argument index to retrieve.
 * @val: pointer to the output value.
 * @sz: output size.
 */
int fn_args_set (fn_arg *argdef, const int i, void *val, size_t sz) {
  /* declare a required variable. */
  int n;

  /* ensure the argdef array index is in bounds. */
  for (n = 0; argdef[n].name; n++) {}
  if (i >= n)
    throw("argument index %d out of bounds [0,%d)", i, n);

  /* determine the type of the argument value. */
  switch (argdef[i].type) {
    /* integer. */
    case FN_VALTYPE_INT:
      argdef[i].val.i = *((int*) val);
      break;

    /* int-array. */
    case FN_VALTYPE_INTS:
      argdef[i].val.iv = *((int**) val);
      break;

    /* boolean. */
    case FN_VALTYPE_BOOL:
      argdef[i].val.b = *((int*) val);
      break;

    /* float. */
    case FN_VALTYPE_FLOAT:
      argdef[i].val.f = *((real*) val);
      break;

    /* float-array. */
    case FN_VALTYPE_FLOATS:
      argdef[i].val.fv = *((real**) val);
      break;

    /* string. */
    case FN_VALTYPE_STRING:
      argdef[i].val.s = *((char**) val);
      break;

    /* chunk. */
    case FN_VALTYPE_CHUNK:
      argdef[i].val.p = *((void**) val);
      break;

    /* other. */
    default:
      throw("unsupported argument type");
  }

  /* store the argument size. */
  argdef[i].sz = sz;

  /* return success. */
  return 1;
}

/* fn_args_get_all(): retrieve the values of every argument in an argument
 * definition array.
 * @argdef: the argument definition array to pull values from.
 * @...: result pointers to store argument values into.
 */
int fn_args_get_all (const fn_arg *argdef, ...) {
  /* declare a few required variables:
   * @vl: variable argument list structure.
   * @val: result pointer of the currently indexed argument.
   * @defi: argument definition array index.
   */
  va_list vl;
  void *val;
  int defi;

  /* initialize the variable arguments list. */
  va_start(vl, argdef);

  /* loop over the argument definition array elements. */
  for (defi = 0; argdef[defi].name; defi++) {
    /* pull another result pointer from the variables arguments list. */
    val = va_arg(vl, void*);

    /* get the value from the array and store it in the result pointer. */
    if (!fn_args_get(argdef, defi, val, NULL))
      throw("failed to get value of argument '%s'", argdef[defi].name);
  }

  /* clean up the variable arguments list. */
  va_end(vl);

  /* return success. */
  return 1;
}

/* fn_args_get_sizes(): retrieve the sizes of every argument in an argument
 * definition array.
 * @argdef: the argument definition array to pull sizes from.
 * @...: result pointers to store argument sizes into.
 */
int fn_args_get_sizes (const fn_arg *argdef, ...) {
  /* declare a few required variables:
   * @vl: variable argument list structure.
   * @sz: size of the currently indexed argument.
   * @defi: argument definition array index.
   */
  va_list vl;
  size_t *sz;
  int defi;

  /* initialize the variable arguments list. */
  va_start(vl, argdef);

  /* loop over the argument definition array elements. */
  for (defi = 0; argdef[defi].name; defi++) {
    /* pull another result pointer from the variable arguments list. */
    sz = va_arg(vl, size_t*);

    /* if the result pointer is non-null, store the currently indexed
     * argument's size into it.
     */
    if (sz)
      *sz = argdef[defi].sz;
  }

  /* clean up the variable arguments list. */
  va_end(vl);

  /* return success. */
  return 1;
}

/* fn_args_parse_int(): parses an integer from a key-value string array.
 * @v: the string array.
 * @c: the string array length.
 */
int fn_args_parse_int (char **v, int c) {
  /* check the array length. */
  if (c != 2)
    throw("r-value required for integer parsing");

  /* parse the value from the string array. */
  return atoi(v[1]);
}

/* fn_args_parse_intarray(): parses an array of integers from a key-value
 * string array.
 * @v: the string array.
 * @c: the string array length.
 * @sz: the size of the result.
 */
int *fn_args_parse_intarray (char **v, int c, size_t *sz) {
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
  strv = strsplit(v[1], ";", &n);

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
  intv = (int*) calloc(n, sizeof(int));
  if (!intv) {
    /* raise an exception and return nothing. */
    raise("failed to allocate int-array");
    return NULL;
  }

  /* loop over the strings in the string array. */
  for (i = 0; i < n; i++)
    intv[i] = atoi(strv[i]);

  /* free the string array */
  strvfree(strv, n);

  /* return the complete array of integers. */
  *sz = n;
  return intv;
}

/* fn_args_parse_bool(): parses a boolean (in the form of an integer) from
 * a key-value string array.
 * @v: the string array.
 * @c: the string array length.
 */
int fn_args_parse_bool (char **v, int c) {
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
    return strbool(v[1]);
  }
  else
    throw("unsupported argument syntax for boolean parsing");

  /* this will never occur. */
  return 0;
}

/* fn_args_parse_float(): parses a floating-point number from a key-value
 * string array.
 * @v: the string array.
 * @c: the string array length.
 */
real fn_args_parse_float (char **v, int c) {
  /* check the array length. */
  if (c != 2)
    throw("r-value required for float parsing");

  /* parse the value from the string array.
   */
  return atof(v[1]);
}

/* fn_args_parse_floatarray(): parses an array of floating-point numbers from
 * a key-value string array.
 * @v: the string array.
 * @c: the string array length.
 * @sz: the size of the result.
 */
real *fn_args_parse_floatarray (char **v, int c, size_t *sz) {
  /* declare a few required variables.
   */
  unsigned int i, n;
  char **strv;
  real *fltv;

  /* check the array length. */
  if (c != 2) {
    /* raise an exception and return nothing. */
    raise("r-value required for float-array parsing");
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
  strv = strsplit(v[1], ";", &n);

  /* check that the split was successful. */
  if (!strv || n < 1) {
    /* raise an exception and return nothing. */
    raise("failed to split float-array string");
    return NULL;
  }

  /* trim whitespace from the string array. */
  strvtrim(strv, n);

  /* allocate an array of floats having the same number of elements as
   * the string array.
   */
  fltv = (real*) calloc(n, sizeof(real));
  if (!fltv) {
    /* raise an exception and return nothing. */
    raise("failed to allocate float-array");
    return NULL;
  }

  /* loop over the strings in the string array. */
  for (i = 0; i < n; i++)
    fltv[i] = atof(strv[i]);

  /* free the string array. */
  strvfree(strv, n);

  /* return the complete array of floats. */
  *sz = n;
  return fltv;
}

/* fn_args_parse_str(): parses a string argument from a key-value string array.
 * this allocates a new string in memory, which needs to be freed by the
 * calling function of fn_scan_args().
 * @v: the string array.
 * @c: the string array length.
 * @sz: the size of the result.
 */
char *fn_args_parse_str (char **v, int c, size_t *sz) {
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
  *sz = n;
  return str;
}

/* fn_args_from_string(): scan a processing function string and extract all
 * known function arguments. when unknown arguments or invalid argument values
 * are supplied, an error is thrown.
 * @argdef: the argument definition array to fill from the string.
 * @argstr: the processing function string to parse.
 */
int fn_args_from_string (fn_arg *argdef, const char *argstr) {
  /* declare variables for argument list parsing:
   * @argv: the argument string array.
   * @argc: the argument string array length.
   * @argi: the argument string array index.
   */
  unsigned int argc, argi;
  char **argv;

  /* declare variables for individual argument parsing:
   * @valv: argument key-value string array.
   * @valc: argument key-value string array length.
   * @defi: argument definition structure loop index.
   */
  unsigned int valc, defi, found;
  char **valv, *noname;
  size_t valsz;

  /* do not attempt to fill a null argdef array. */
  if (!argdef)
    return 1;

  /* parse the argument string into an array. */
  argv = strsplit(argstr, ",", &argc);

  /* check that the split was successful. */
  if (!argv || !argc)
    throw("failed to split argument string");

  /* trim the string array elements. */
  strvtrim(argv, argc);

  /* loop over the argument definition entries. */
  for (defi = 0; argdef[defi].name; defi++) {
    /* initialize the negated argument name. */
    noname = NULL;

    /* check if a negated argument should be searched for as well. */
    if (argdef[defi].type == FN_VALTYPE_BOOL) {
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

    /* check if a value was found for the current argument. */
    if (found) {
      /* initialize the argument value size. */
      valsz = 0;

      /* act based on the expected type of the argument. */
      switch (argdef[defi].type) {
        /* signed int. */
        case FN_VALTYPE_INT:
          argdef[defi].val.i = fn_args_parse_int(valv, valc);
          break;

        /* int array. */
        case FN_VALTYPE_INTS:
          argdef[defi].val.iv = fn_args_parse_intarray(valv, valc, &valsz);
          break;

        /* boolean. */
        case FN_VALTYPE_BOOL:
          argdef[defi].val.b = fn_args_parse_bool(valv, valc);
          break;

        /* float. */
        case FN_VALTYPE_FLOAT:
          argdef[defi].val.f = fn_args_parse_float(valv, valc);
          break;

        /* float array. */
        case FN_VALTYPE_FLOATS:
          argdef[defi].val.fv = fn_args_parse_floatarray(valv, valc, &valsz);
          break;

        /* string. */
        case FN_VALTYPE_STRING:
          argdef[defi].val.s = fn_args_parse_str(valv, valc, &valsz);
          break;

        /* other... */
        default:
          throw("unsupported argument type '%c'", argdef[defi].type);
      }

      /* free the values string array. */
      strvfree(valv, valc);

      /* store the argument value size. */
      argdef[defi].sz = valsz;
    }

    /* free the negated value string, if required. */
    if (noname)
      free(noname);
  }

  /* free the argument string array. */
  strvfree(argv, argc);

  /* return success. */
  return 1;
}

/* fn_args_copy(): copy an argument definition (argdef) array to a newly
 * allocated location in memory.
 * @argsrc: source argdef array to duplicate.
 */
fn_arg *fn_args_copy (const fn_arg *argsrc) {
  /* declare a few required variables. */
  fn_arg *args;
  int n;

  /* return nothing if the source argdef array is null. */
  if (!argsrc)
    return NULL;

  /* compute the size of the argdef array. */
  for (n = 0; argsrc[n].name; n++) {}
  n = (n + 1) * sizeof(fn_arg);

  /* allocate memory for the array duplicate. */
  args = (fn_arg*) malloc(n);

  /* ensure allocation was successful. */
  if (!args) {
    /* raise an error and return nothing. */
    raise("failed to allocate duplicate argdef array");
    return NULL;
  }

  /* copy the argdef contents onto the allocated memory. */
  memcpy(args, argsrc, n);

  /* return the newly allocated, duplicated argdef array. */
  return args;
}

