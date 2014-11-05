
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

/* include the nmr datum header. */
#include <hxnd/nmr-datum.h>

/* define argument type characters used when parsing processing function
 * strings.
 */
#define FN_ARGTYPE_INT    'i'
#define FN_ARGTYPE_BOOL   'b'
#define FN_ARGTYPE_FLOAT  'f'

/* define string names for all available processing functions.
 */
#define FN_NAME_FFT  "fft"
#define FN_NAME_HT   "ht"

/* fn_args: structure for holding information on arguments accepted by any
 * given processing function.
 */
typedef struct {
  /* @name: argument name string.
   * @type: argument type character.
   */
  const char *name;
  char type;

  /* @def: default argument value (string-representation).
   */
  const char *def;
}
fn_args;

/* function declarations: */

int fn_scan_args (const char *argstr, const fn_args *argdef, ...);

/* function declarations, processing: */

int fn_execute (datum *D, const char *name, const int dim, const char *argstr);

int fn_execute_fft (datum *D, const int dim, const char *argstr);

int fn_execute_ht (datum *D, const int dim, const char *argstr);

#endif /* __HXND_FN_H__ */

