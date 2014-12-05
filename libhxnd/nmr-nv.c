
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
  /* declare a few required variables:
   * @fh: input file handle.
   * @wd: first word of the file.
   */
  int32_t wd;
  FILE *fh;

  /* open the input file. */
  fh = fopen(fname, "rb");

  /* check that the file was opened. */
  if (!fh)
    throw("failed to open '%s'", fname);

  /* read the first word. */
  if (!fread(&wd, sizeof(int32_t), 1, fh))
    throw("failed to read magic number");

  /* close the input file. */
  fclose(fh);

  /* check the magic word, without swapping. */
  if (wd == NV_MAGIC)
    return 1;

  /* swap the word and check again. */
  bytes_swap_u32((uint32_t*) &wd);
  if (wd == NV_MAGIC)
    return 1;

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
  /* declare a few required variables:
   */
  unsigned int i, n_bytes, n_words;
  uint8_t *bytes;

  /* read in the file header bytes. */
  n_bytes = sizeof(struct nv_file_header);
  n_words = n_bytes / sizeof(int32_t);
  bytes = bytes_read_block(fname, 0, n_bytes);

  /* check that the bytes were read successfully. */
  if (!bytes)
    throw("failed to read header from '%s'", fname);

  /* copy the header bytes onto the header structure. */
  memcpy(hdr, bytes, n_bytes);

  /* check if the magic word is correct. if not, the bytes likely need
   * to be swapped.
   */
  if (hdr->magic != NV_MAGIC) {
    /* swap all header bytes. */
    bytes_swap(bytes, n_words, sizeof(int32_t));

    /* re-copy the bytes onto the structure. */
    memcpy(hdr, bytes, n_bytes);

    /* swap the file header strings back. */
    bytes_swap((uint8_t*) hdr->sequence, NV_HDRSTR_SZ_SEQUENCE, sizeof(char));
    bytes_swap((uint8_t*) hdr->comment, NV_HDRSTR_SZ_COMMENT, sizeof(char));

    /* swap the dimension header strings back. */
    for (i = 0; i < NV_MAXDIM; i++)
      bytes_swap((uint8_t*) hdr->dims[i].label,
                 NV_HDRSTR_SZ_LABEL, sizeof(char));

    /* opposite endianness. */
    *endianness = bytes_get_nonnative();
  }
  else {
    /* same endianness. */
    *endianness = bytes_get_native();
  }

  /* free the read bytes. */
  free(bytes);

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
  /* declare variables required to determine byte ordering:
   * @endianness: the byte ordering of the data file.
   * @hdr: the nmrview file header structure.
   */
  enum byteorder endianness = BYTES_ENDIAN_AUTO;
  struct nv_file_header hdr;

  /* read the header information from the data file. */
  if (!nv_read_header(fname, &endianness, &hdr))
    throw("failed to read header of '%s'", fname);

  /* FIXME: implement nv_fill_datum() */

  /* store the filename string. */
  D->fname = (char*) malloc((strlen(fname) + 1) * sizeof(char));
  if (D->fname)
    strcpy(D->fname, fname);

  /* store the datum type. */
  D->type = DATUM_TYPE_NV;
  D->endian = endianness;

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

