
/* ndmath: A framework for loading free induction decay series data.
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
#ifndef __NDMATH_SER_H__
#define __NDMATH_SER_H__

/* include the n-dimensional math header. */
#include <hxnd/hx.h>

/* function declarations: */

int ser_toarray (uint8_t *bytes, unsigned int nbytes,
                 int endianness, int wordsz, int flt,
                 hx_array *x);

int ser_deinterlace (hx_array *x);

#endif /* __NDMATH_SER_H__ */

