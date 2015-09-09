
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
#ifndef __GHX_GHX_APP_WINDOW_H__
#define __GHX_GHX_APP_WINDOW_H__

/* include the gtk library header. */
#include <gtk/gtk.h>

/* define a macro that returns the GhxAppWindow GType. */
#define GHX_APP_WINDOW_TYPE (ghx_app_window_get_type())

/* define a macro that casts a GObject as a GhxAppWindow. */
#define GHX_APP_WINDOW(obj) \
  (G_TYPE_CHECK_INSTANCE_CAST((obj), GHX_APP_WINDOW_TYPE, GhxAppWindow))

/* declare the required structures. */
typedef struct _GhxAppWindow        GhxAppWindow;
typedef struct _GhxAppWindowClass   GhxAppWindowClass;
typedef struct _GhxAppWindowPrivate GhxAppWindowPrivate;

/* function declarations: */

GType ghx_app_window_get_type (void);

GhxAppWindow *ghx_app_window_new (GhxApp *app);

void ghx_app_window_open (GhxAppWindow *win, GFile *file);

#endif /* __GHX_GHX_APP_WINDOW_H__ */

