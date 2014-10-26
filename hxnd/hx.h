
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
#ifndef __HXND_HX_H__
#define __HXND_HX_H__

/* include required standard c library headers. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

/* define 'pi' if needed. */
#ifndef M_PI
#define M_PI 3.1415926539
#endif

/* real: type definition of a real scalar floating-point value.
 */
typedef float real;

/* include all hypercomplex math headers. */
#include "hx-algebra.h"
#include "hx-scalar.h"
#include "hx-index.h"
#include "hx-array.h"
#include "hx-cmp.h"
#include "hx-arith.h"
#include "hx-fourier.h"

#endif /* __HXND_HX_H__ */

