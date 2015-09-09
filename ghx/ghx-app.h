
/* hxnd: A framework for n-dimensional hypercomplex calculations for NMR.
 * Copyright (C) 2014-2015  Bradley Worley  <geekysuavo@gmail.com>.
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
#ifndef __GHX_GHX_APP_H__
#define __GHX_GHX_APP_H__

/* include the gtk library header. */
#include <gtk/gtk.h>

/* define a macro that returns the GhxApp GType. */
#define GHX_APP_TYPE (ghx_app_get_type())

/* define a macro that casts a GObject as a GhxApp. */
#define GHX_APP(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST((obj), GHX_APP_TYPE, GhxApp))

/* declare the required structures. */
typedef struct _GhxApp      GhxApp;
typedef struct _GhxAppClass GhxAppClass;

/* function declarations: */

GType ghx_app_get_type (void);

GhxApp *ghx_app_new (void);

#endif /* __GHX_GHX_APP_H__ */

