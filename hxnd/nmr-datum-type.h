
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
#ifndef __HXND_NMR_DATUM_TYPE_H__
#define __HXND_NMR_DATUM_TYPE_H__

/* include the nmr format headers. */
#include <hxnd/nmr-hxnd.h>
#include <hxnd/nmr-text.h>
#include <hxnd/nmr-bruker.h>
#include <hxnd/nmr-varian.h>
#include <hxnd/nmr-pipe.h>
#include <hxnd/nmr-ucsf.h>
#include <hxnd/nmr-nv.h>
#include <hxnd/nmr-rnmrtk.h>

/* function declarations (nmr-datum-type.c): */

const char *datum_type_name (enum datum_type typ);

const char *datum_type_desc (enum datum_type typ);

enum datum_type datum_type_lookup (const char *name);

enum datum_type datum_type_guess (const char *fname);

int datum_type_decode (datum *D, const char *fname, enum datum_type typ);

int datum_type_encode (datum *D, const char *fname, enum datum_type typ);

int datum_type_array (datum *D, enum datum_type typ);

int datum_type_post (datum *D, enum datum_type typ);

#endif /* __HXND_NMR_DATUM_TYPE_H__ */

