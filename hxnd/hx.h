
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
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* include the openmp library header. */
#include <omp.h>

/* define 'pi' if needed. */
#ifndef M_PI
#define M_PI 3.1415926539
#endif

/* include the definition of real scalar values. */
#include <hxnd/hx-real.h>

/* include the traceback header. */
#include <hxnd/trace.h>

/* include all hypercomplex math headers. */
#include <hxnd/hx-algebra.h>
#include <hxnd/hx-scalar.h>
#include <hxnd/hx-index.h>
#include <hxnd/hx-array.h>
#include <hxnd/hx-cmp.h>
#include <hxnd/hx-arith.h>
#include <hxnd/hx-blas.h>
#include <hxnd/hx-phasor.h>
#include <hxnd/hx-fourier.h>
#include <hxnd/hx-baseline.h>
#include <hxnd/hx-ist.h>
#include <hxnd/hx-maxent.h>

#endif /* __HXND_HX_H__ */

