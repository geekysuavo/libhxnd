
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
#ifndef __HXND_TRACE_H__
#define __HXND_TRACE_H__

/* include required standard c headers. */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>

/* raise(): macro function to add traceback information onto a stack trace.
 */
#define raise(...) \
  traceback_throw(__FILE__, __LINE__, __VA_ARGS__)

/* throw(): just like raise(), but this macro returns the result of
 * traceback_throw(), which will be zero.
 */
#define throw(...) \
  return traceback_throw(__FILE__, __LINE__, __VA_ARGS__)

/* trace(): just like raise, but this macro also prints the resulting stack
 * trace and returns a value that indicates an error in main() functions.
 */
#define trace(...) \
  { traceback_throw(__FILE__, __LINE__, __VA_ARGS__); \
    traceback_print(); \
    return 1; }

/* function declarations: */

void traceback_init (void);

void traceback_print (void);

void traceback_clear (void);

int traceback_throw (const char *f, const unsigned int l,
                     const char *format, ...);

#endif /* __HXND_TRACE_H__ */

