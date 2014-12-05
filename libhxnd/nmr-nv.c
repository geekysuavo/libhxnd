
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

/* include the nmrview header. */
#include <hxnd/nmr-nv.h>

/* nv_check_magic(): checks the first "magic" bytes of a file and
 * returns whether they match the nmrview-format file magic number.
 * @fname: the input filename.
 */
int nv_check_magic (const char *fname) {
  /* FIXME: implement nv_check_magic() */

  /* no match. */
  return 0;
}

/* nv_read_header(): read the contents of a header from an nmrview-format file.
 * @fname: the input filename.
 * @endianness: byte order result pointer.
 * @hdr: file header result pointer.
 */
int nv_read_header (const char *fname,
                    enum byteorder *endianness,
                    struct nv_file_header *hdr) {
  /* FIXME: implement nv_read_header() */

  /* return success. */
  return 1;
}

/* nv_read(): reads an nmrview-format data file into a real linear array.
 * @fname: the input data filename.
 * @x: the output array.
 */
int nv_read (const char *fname, hx_array *x) {
  /* FIXME: implement nv_read() */

  /* return success. */
  return 1;
}

/* nv_fill_datum(): intelligently parses nmrview parameters into an
 * NMR datum structure.
 * @fname: the input filename.
 * @D: pointer to the datum struct to fill.
 */
int nv_fill_datum (const char *fname, datum *D) {
  /* FIXME: implement nv_fill_datum() */

  /* return success. */
  return 1;
}

/* nv_fwrite_datum(): writes an NMR datum structure in nmrview-format to
 * an opened file stream.
 * @D: pointer to the source structure.
 * @fh: the output file stream.
 */
int nv_fwrite_datum (datum *D, FILE *fh) {
  /* FIXME: implement nv_fwrite_datum() */

  /* return success. */
  return 1;
}

