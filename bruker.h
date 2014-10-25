
/* ndmath: A framework for n-dimensional hypercomplex calculations for NMR
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
#ifndef __HXND_BRUKER_H__
#define __HXND_BRUKER_H__

/* include the n-dimensional math header. */
#include <hxnd/hx.h>

/* include the byte loading and serial headers. */
#include "bytes.h"
#include "ser.h"

/* define constant parameter type characters for acqus file parsing.
 */
#define BRUKER_PARMTYPE_INT     'i'
#define BRUKER_PARMTYPE_FLOAT   'f'
#define BRUKER_PARMTYPE_STRING  's'

/* function declarations: */

int bruker_read_parms (const char *fname, unsigned int n, ...);

int bruker_read (const char *fname, unsigned int endianness,
                 unsigned int nblk, unsigned int szblk,
                 hx_array *x);

#endif /* __HXND_BRUKER_H__ */

