
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
#ifndef __HXND_STR_H__
#define __HXND_STR_H__

/* include required standard c headers. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* function declarations: */

void strtrim (char *s);

void strltrim (char *s);

void strrtrim (char *s);

void strnltrim (char *s);

int strbool (const char *s);

char **strsplit (const char *s1, const char *s2, unsigned int *ntok);

void strvtrim (char **strv, unsigned int n);

void strvfree (char **strv, unsigned int n);

#endif /* __HXND_STR_H__ */

