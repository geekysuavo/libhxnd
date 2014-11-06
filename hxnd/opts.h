
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
#ifndef __HXND_OPTS_H__
#define __HXND_OPTS_H__

/* include the string handling header. */
#include <hxnd/str.h>

/* opts_def: options definition structure for informing opts_get() about
 * all possible command-line options that are accepted by a program.
 */
typedef struct {
  /* @lname: long option name string.
   * @sname: short option name char.
   * @has_arg: whether this option requires an argument.
   */
  const char *lname;
  int has_arg;
  char sname;
}
opts_def;

/* function declarations: */

int opts_get (int argc, char **argv, const opts_def *opts, int *argi);

int opts_parse_arg (char *arg, const char *delim,
                    char **lvalue, char **rvalue,
                    int *d);

#endif /* __HXND_OPTS_H__ */

