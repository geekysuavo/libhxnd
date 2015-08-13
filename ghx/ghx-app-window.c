
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

/* include the graphical interface header. */
#include <ghx/ghx.h>

#define GHX_APP_WINDOW_GET_PRIVATE(obj) \
  (G_TYPE_INSTANCE_GET_PRIVATE((obj), GHX_APP_WINDOW_TYPE, \
                               GhxAppWindowPrivate))

/* define the GhxAppWindow object struct. */
struct _GhxAppWindow {
  /* parent window object. */
  GtkApplicationWindow parent;

  /* private data structure. */
  GhxAppWindowPrivate *priv;
};

/* define the GhxAppWindow class struct. */
struct _GhxAppWindowClass {
  /* parent window class. */
  GtkApplicationWindowClass parent_class;
};

/* define the GhxAppWindow private struct. */
struct _GhxAppWindowPrivate {
  /* the main header bar. */
  GtkHeaderBar *header;
};

/* define the GhxAppWindow and its core functions. */
G_DEFINE_TYPE(GhxAppWindow, ghx_app_window, GTK_TYPE_APPLICATION_WINDOW);

/* ghx_app_init_window(): initialize a GhxAppWindow instance.
 * @win: GhxAppWindow instance to initialize.
 */
static void ghx_app_window_init (GhxAppWindow *win) {
  /* initialize the window's private data structure. */
  win->priv = GHX_APP_WINDOW_GET_PRIVATE(win);

  /* create a header bar. */
  win->priv->header = GTK_HEADER_BAR(gtk_header_bar_new());

  /* set the initial header bar title and subtitle. */
  gtk_header_bar_set_title(win->priv->header, "ghx");
  gtk_header_bar_set_subtitle(win->priv->header, "~");

  /* show the shrink/expand/close button group on the header bar. */
  gtk_header_bar_set_show_close_button(win->priv->header, TRUE);

  /* attach and show the header bar. */
  gtk_widget_show(GTK_WIDGET(win->priv->header));
  gtk_window_set_titlebar(GTK_WINDOW(win),
                          GTK_WIDGET(win->priv->header));

  /* FIXME: implement window initialization here. */

  /* set the default size of the window. */
  gtk_window_set_default_size(GTK_WINDOW(win), 600, 400);
}

/* ghx_app_window_class_init(): initialize a GhxAppWindowClass instance.
 * @class: GhxAppWindowClass instance to initialize.
 */
static void ghx_app_window_class_init (GhxAppWindowClass *class) {
  /* register the private data structure with the class. */
  g_type_class_add_private(G_OBJECT_CLASS(class),
                           sizeof(GhxAppWindowPrivate));
}

/* ghx_app_window_new(): create a new GhxAppWindow instance.
 * @app: GhxApp instance to assign the window to.
 */
GhxAppWindow *ghx_app_window_new (GhxApp *app) {
  /* create and return the new window object. */
  return g_object_new(GHX_APP_WINDOW_TYPE,
                      "application", app,
                      "show-menubar", FALSE,
                      "type", GTK_WINDOW_TOPLEVEL,
                      NULL);
}

/* ghx_app_window_open(): open a file with a GhxAppWindow instance.
 * @win: GhxAppWindow instance to utilize.
 * @file: GFile to open with @win.
 */
void ghx_app_window_open (GhxAppWindow *win, GFile *file) {
  /* FIXME: implement ghx_app_window_open() */
}

