
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

/* include the graphical interface header. */
#include <ghx/ghx.h>

/* define the GhxApp object struct. */
struct _GhxApp {
  /* parent application object. */
  GtkApplication parent;
};

/* define the GhxApp class struct. */
struct _GhxAppClass {
  /* parent application class. */
  GtkApplicationClass parent_class;
};

/* define the GhxApp and its core functions. */
G_DEFINE_TYPE(GhxApp, ghx_app, GTK_TYPE_APPLICATION);

/* ghx_app_init(): initialize a GhxApp instance.
 * @app: GhxApp instance to initialize.
 */
static void ghx_app_init (GhxApp *app) {
  /* initialization tasks will be placed here. */
}

/* ghx_app_activate(): activate a GhxApp without input files.
 * @app: application instance to activate.
 */
static void ghx_app_activate (GApplication *app) {
  /* declare required variables:
   * @win: primary application window.
   */
  GhxAppWindow *win;

  /* create and present new window. */
  win = ghx_app_window_new(GHX_APP(app));
  gtk_window_present(GTK_WINDOW(win));
}

/* ghx_app_open(): activate a GhxApp with a list of input files.
 * @app: application instance to activate.
 * @files: array of GFile pointers.
 * @n_files: length of the @files array.
 * @hint: hint provided by the instance.
 */
static void ghx_app_open (GApplication *app,
                          GFile **files,
                          gint n_files,
                          const gchar *hint) {
  /* declare required variables:
   * @win: primary application window.
   * @winlst: list of open application windows.
   * @i: loop counter and file index.
   */
  GhxAppWindow *win;
  GList *winlst;
  int i;

  /* get the list of windows owned by the application instance. */
  winlst = gtk_application_get_windows(GTK_APPLICATION(app));

  /* check if the list of windows is empty. */
  if (winlst) {
    /* non-empty list. access the window. */
    win = GHX_APP_WINDOW(winlst->data);
  }
  else {
    /* empty list. create a new window. */
    win = ghx_app_window_new(GHX_APP(app));
  }

  /* loop over each file to be opened. */
  for (i = 0; i < n_files; i++) {
    /* open the currently indexed file. */
    ghx_app_window_open(win, files[i]);
  }

  /* present the primary application window. */
  gtk_window_present(GTK_WINDOW(win));
}

/* ghx_app_class_init(): initialize a GhxAppClass instance.
 * @class: GhxAppClass instance to initialize.
 */
static void ghx_app_class_init (GhxAppClass *class) {
  /* link the "activate" and "open" function pointer into the class. */
  G_APPLICATION_CLASS(class)->activate = ghx_app_activate;
  G_APPLICATION_CLASS(class)->open = ghx_app_open;
}

/* ghx_app_new(): create a new GhxApp instance.
 */
GhxApp *ghx_app_new (void) {
  /* create and return the new application object. */
  return g_object_new(GHX_APP_TYPE,
                      "application-id", "org.gtk.ghx",
                      "flags", G_APPLICATION_HANDLES_OPEN,
                      NULL);
}

