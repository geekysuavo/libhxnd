
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
#ifndef __HXND_FN_HANDLERS_H__
#define __HXND_FN_HANDLERS_H__

/* include the nmr datum and dataset headers. */
#include <hxnd/nmr-datum.h>
#include <hxnd/mx-dataset.h>

/* define string names for all available processing functions.
 */
#define FN_NAME_ABS       "abs"
#define FN_NAME_ADD       "add"
#define FN_NAME_COMPLEX   "complex"
#define FN_NAME_CUT       "cut"
#define FN_NAME_FFT       "fft"
#define FN_NAME_HT        "ht"
#define FN_NAME_IST       "ist"
#define FN_NAME_MIRROR    "mirror"
#define FN_NAME_MULTIPLY  "multiply"
#define FN_NAME_PHASE     "phase"
#define FN_NAME_PROJECT   "project"
#define FN_NAME_REAL      "real"
#define FN_NAME_RESIZE    "resize"
#define FN_NAME_SHIFT     "shift"
#define FN_NAME_WINDOW    "window"
#define FN_NAME_ZEROFILL  "zerofill"

/* define string names for all available multivariate functions.
 */
#define FN_NAME_NORM   "norm"
#define FN_NAME_SCALE  "scale"

/* function declarations, datum handlers: */

int fn_abs (datum *D, const int dim, const fn_arg *args);

int fn_add (datum *D, const int dim, const fn_arg *args);

int fn_complex (datum *D, const int dim, const fn_arg *args);

int fn_cut (datum *D, const int dim, const fn_arg *args);

int fn_fft (datum *D, const int dim, const fn_arg *args);

int fn_ht (datum *D, const int dim, const fn_arg *args);

int fn_ist (datum *D, const int dim, const fn_arg *args);

int fn_mirror (datum *D, const int dim, const fn_arg *args);

int fn_multiply (datum *D, const int dim, const fn_arg *args);

int fn_phase (datum *D, const int dim, const fn_arg *args);

int fn_project (datum *D, const int dim, const fn_arg *args);

int fn_real (datum *D, const int dim, const fn_arg *args);

int fn_resize (datum *D, const int dim, const fn_arg *args);

int fn_shift (datum *D, const int dim, const fn_arg *args);

int fn_window (datum *D, const int dim, const fn_arg *args);

int fn_zerofill (datum *D, const int dim, const fn_arg *args);

/* function declarations, dataset handlers: */

int fn_norm (dataset *Dset, const int dim, const fn_arg *args);

int fn_scale (dataset *Dset, const int dim, const fn_arg *args);

#endif /* __HXND_FN_HANDLERS_H__ */

