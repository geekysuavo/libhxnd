
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

/* include the rnmrtk header. */
#include <hxnd/nmr-rnmrtk.h>

/* declare a buffer size for reading rnmrtk parameter files.
 */
#define N_BUF  256

/* rnmrtk_check_file(): check a file for the requisite properties that may
 * define it as an RNMRTK data file.
 * @fname: the input data file name.
 */
int rnmrtk_check_file (const char *fname) {
  /* FIXME: implement rnmrtk_check_file() */
  return 0;
}

/* rnmrtk_read_parms(): read the parameter file corresponding to an RNMRTK
 * data file.
 * @fname: the data filename.
 * @par: pointer to the output parameter structure.
 */
int rnmrtk_read_parms (const char *fname, struct rnmrtk_parms *par) {
  /* FIXME: implement rnmrtk_read_parms() */
  throw("unimplemented!");
}

/* rnmrtk_read(): read an RNMRTK data file into a real linear array.
 * @fname: the input data filename.
 * @x: the output array.
 */
int rnmrtk_read (const char *fname, hx_array *x) {
  /* FIXME: implement rnmrtk_read() */
  throw("unimplemented!");
}

/* rnmrtk_fill_datum(): load parameter information accompanying an RNMRTK
 * data file into an NMR datum structure.
 * @fname: the input data file name.
 * @D: pointer to the datum struct to fill.
 */
int rnmrtk_fill_datum (const char *fname, datum *D) {
  /* FIXME: implement rnmrtk_fill_datum() */
  throw("unimplemented!");
}

