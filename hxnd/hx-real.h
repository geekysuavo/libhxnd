
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
#ifndef __HXND_HX_REAL_H__
#define __HXND_HX_REAL_H__

/* HX_DOUBLE_PRECISION: preprocessor declaration indicating that all real
 * scalar floating-point values are to be double-precision.
 */
// #define HX_DOUBLE_PRECISION  1

/* real: type definition of a real scalar floating-point value.
 */
#ifdef HX_DOUBLE_PRECISION
  /* use double-precision floats. */
  typedef double real;
#else
  /* use single-precision floats. */
  typedef float real;
#endif

#endif /* __HXND_HX_REAL_H__ */

