
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

/* include the pipe header. */
#include <hxnd/nmr-pipe.h>

/* pipe_read_header(): read the contents of a header from a pipe-format file.
 * @fname: the input filename.
 * @endianness: byte order result pointer.
 * @hdr: header result pointer.
 */
int pipe_read_header (const char *fname, enum byteorder *endianness,
                      struct pipe_header *hdr) {
  /* declare a few required variables:
   */
  unsigned int n_bytes;
  uint8_t *bytes;

  /* read in the file header bytes. */
  n_bytes = sizeof(struct pipe_header);
  bytes = bytes_read_block(fname, 0, n_bytes);

  /* check that the bytes were read successfully. */
  if (!bytes)
    throw("failed to read header from '%s'", fname);

  /* copy the header bytes onto the header structure. */
  memcpy(hdr, bytes, n_bytes);

  /* check if the floating point order value is correct. if not, the bytes
   * likely need to be swapped.
   */
  if (hdr->order != 0.0 && hdr->order != (float) 2.345) {
    /* swap the bytes of each word. */
    bytes_swap_general(bytes, n_bytes, sizeof(float));

    /* re-copy the bytes onto the structure. */
    memcpy(hdr, bytes, n_bytes);

    /* opposite endianness. */
    if (bytes_native(BYTES_ENDIAN_BIG))
      *endianness = BYTES_ENDIAN_LITTLE;
    else
      *endianness = BYTES_ENDIAN_BIG;
  }
  else {
    /* same endianness. */
    if (bytes_native(BYTES_ENDIAN_BIG))
      *endianness = BYTES_ENDIAN_BIG;
    else
      *endianness = BYTES_ENDIAN_LITTLE;
  }

  /* free the read bytes. */
  free(bytes);

  /* return success. */
  return 1;
}

/* pipe_fill_datum(): intelligently parses pipe parameters into an
 * NMR datum structure.
 * @fname: the input filename.
 * @D: pointer to the datum struct to fill.
 */
int pipe_fill_datum (const char *fname, datum *D) {
  /* declare variables required to determine byte ordering. */
  enum byteorder endianness = BYTES_ENDIAN_AUTO;
  struct pipe_header hdr;

  /* read the header information from the data file. */
  if (!pipe_read_header(fname, &endianness, &hdr))
    throw("failed to read header of '%s'", fname);

  /* FIXME: implement pipe_fill_datum() */

  /* store the filename string. */
  D->fname = (char*) malloc((strlen(fname) + 1) * sizeof(char));
  if (D->fname)
    strcpy(D->fname, fname);

  /* store the datum type. */
  D->type = DATUM_TYPE_PIPE;
  D->endian = endianness;

  /* return success. */
  return 1;
}

