
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

/* include the option parsing and traceback headers. */
#include <hxnd/opts.h>
#include <hxnd/trace.h>

/* opts_get(): returns the next parsed option in an argument array.
 * this function does not implement any kind of compliant option parsing,
 * but instead just allows for single-character options to be parsed,
 * without the ability to group options into one argument. Long options
 * are also supported.
 * @argc: argument count.
 * @argv: argument array.
 * @opts: long option definition struct array.
 * @argi: pointer to current argument index.
 */
int opts_get (int argc, char **argv, const opts_def *opts, int *argi) {
  /* declare a few required variables. */
  char *arg;
  int i;

  /* ensure parsing begins at the correct index. */
  if (*argi < 1)
    *argi = 1;

  /* check that the index is in bounds. */
  if (*argi >= argc)
    return -1;

  /* get the currently indexed option. */
  arg = argv[*argi];

  /* check the currently indexed option. */
  if (strlen(arg) == 2 && arg[0] == '-') {
    /* option char identified. search the option definition array. */
    for (i = 0; opts[i].lname; i++) {
      /* check if the current option is a match. */
      if (opts[i].sname == arg[1])
        break;
    }

    /* check if an option was identified. */
    if (opts[i].lname) {
      /* yes. check if an option argument is required. */
      if (opts[i].has_arg) {
        /* ensure an extra argument exists. */
        if (*argi >= argc - 1)
          throw("option '-%c' requires an argument", arg[1]);

        /* increment the argument index to skip the option argument. */
        (*argi)++;
      }

      /* everything checks out. increment the argument index and return. */
      (*argi)++;
      return (int) arg[1];
    }
    else {
      /* throw an exception. */
      throw("invalid short option '-%c'", arg[1]);
    }
  }
  else if (strlen(arg) >= 3 && arg[0] == '-' && arg[1] == '-') {
    /* option string identified. search the option definition array. */
    for (i = 0; opts[i].lname; i++) {
      /* check if the current option is a match. */
      if (strcmp(opts[i].lname, arg + 2) == 0)
        break;
    }

    /* check if an option was identified. */
    if (opts[i].lname) {
      /* yes. check if an option argument is required. */
      if (opts[i].has_arg) {
        /* ensure an extra argument exists. */
        if (*argi >= argc - 1)
          throw("option '--%s' requires an argument", arg + 2);

        /* increment the argument index to skip the option argument. */
        (*argi)++;
      }

      /* everything checks out. increment the argument index and return. */
      (*argi)++;
      return (int) opts[i].sname;
    }
    else {
      /* throw an exception. */
      throw("invalid long option '--%s'", arg + 2);
    }
  }
  else {
    /* throw an exception. */
    throw("invalid argument '%s'", arg);
  }

  /* return nothing found. */
  return -1;
}

