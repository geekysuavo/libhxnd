
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

/* nv_tiler(): redo or undo the tiling inherent in nmrview-format files,
 * effectively mapping between tiled array data and a real linear
 * array suitable for datum_refactor_array().
 * @x: pointer to the array to linearize.
 * @hdr: pointer to the file header.
 * @dir: direction, either 0 (linearize) or 1 (tileize).
 */
int nv_tiler (hx_array *x, struct nv_file_header *hdr, int dir) {
  /* declare a few required variables:
   * @i: general purpose loop counter.
   * @k: number of dimensions.
   * @nt: array of tile counts.
   * @szt: array of tile sizes.
   */
  int i, k, *nt, *szt;

  /* gain a handle on the dimensionality of the array data. */
  k = (int) hdr->ndims;

  /* allocate the size index arrays. */
  nt = hx_array_index_alloc(k);
  szt = hx_array_index_alloc(k);

  /* check that all index allocations were successful. */
  if (!nt || !szt)
    throw("failed to allocate two sets of %d indices", k);

  /* initialize the size arrays. */
  for (i = 0; i < k; i++) {
    /* check that the point and tile counts are nonzero. */
    if (hdr->dims[i].sz < 1 || hdr->dims[i].szblk < 1)
      throw("invalid tiling (%d, %d) along dimension %d",
            hdr->dims[i].szblk, hdr->dims[i].sz, i);

    /* check that the point count is evenly divided by the tile size. */
    if (hdr->dims[i].sz % hdr->dims[i].szblk)
      throw("tile size %u does not evenly divide point count %u",
            hdr->dims[i].szblk, hdr->dims[i].sz);

    /* set the tile count. */
    nt[i] = hdr->dims[i].sz / hdr->dims[i].szblk;

    /* set the tile point count. */
    szt[i] = hdr->dims[i].szblk;
  }

  /* perform the mapping operation. */
  if (!hx_array_tiler(x, k, nt, szt, dir))
    throw("failed to perform tile mapping");

  /* free the allocated index arrays. */
  free(nt);
  free(szt);

  /* return success. */
  return 1;
}

/* nv_linearize(): maps tiles to linear array data.
 */
#define nv_linearize(x, hdr) \
  nv_tiler(x, hdr, HX_ARRAY_TILER_FORWARD)

/* ucsf_tileize(): maps linear array data to tiles.
 */
#define nv_tileize(x, hdr) \
  nv_tiler(x, hdr, HX_ARRAY_TILER_REVERSE)

/* nv_read(): reads an nmrview-format data file into a real linear array.
 * @fname: the input data filename.
 * @x: the output array.
 */
int nv_read (const char *fname, hx_array *x) {
  /* declare variables required to read header information:
   * @endianness: the byte ordering of the data file.
   * @hdr: the nmrview file header structure.
   */
  enum byteorder endian = BYTES_ENDIAN_AUTO;
  struct nv_file_header hdr;

  /* declare variables for reading raw data bytes:
   * @offset: byte offset where point data begins.
   * @i: general purpose loop counter.
   * @n: number of words to read.
   */
  unsigned int offset, i, n;
  FILE *fh;

  /* read the header information from the data file. */
  if (!nv_read_header(fname, &endian, &hdr))
    throw("failed to read header of '%s'", fname);

  /* compute the byte offset from which to begin reading point data. */
  offset = sizeof(struct nv_file_header);

  /* compute the number of words to read. */
  for (i = 0, n = 1; i < hdr.ndims; i++)
    n *= hdr.dims[i].sz;

  /* open the input file for reading. */
  fh = fopen(fname, "rb");

  /* check that the file was opened. */
  if (!fh)
    throw("failed to open '%s'", fname);

  /* read data from the file into the output array. */
  if (!hx_array_fread_raw(fh, x, endian, sizeof(float), 1,
                          offset, 0, 1, n, 0))
    throw("failed to read raw data from '%s'", fname);

  /* close the input file. */
  fclose(fh);

  /* convert tiles to linear values. */
  if (!nv_linearize(x, &hdr))
    throw("failed to linearize tiled array");

  /* return success. */
  return 1;
}

/* nv_fill_datum(): intelligently parses nmrview parameters into an
 * NMR datum structure.
 * @fname: the input filename.
 * @D: pointer to the datum struct to fill.
 */
int nv_fill_datum (const char *fname, datum *D) {
  /* declare variables required to read header information:
   * @endianness: the byte ordering of the data file.
   * @hdr: the nmrview file header structure.
   * @d: dimension loop counter.
   */
  enum byteorder endianness = BYTES_ENDIAN_AUTO;
  struct nv_file_header hdr;
  unsigned int d;

  /* read the header information from the data file. */
  if (!nv_read_header(fname, &endianness, &hdr))
    throw("failed to read header of '%s'", fname);

  /* store the dimension count. */
  D->nd = (unsigned int) hdr.ndims;

  /* check the dimension count. */
  if (D->nd < 1 || D->nd > NV_MAXDIM)
    throw("invalid dimensionality %d", hdr.ndims);

  /* allocate the dimension parameter array. */
  D->dims = (datum_dim*) calloc(D->nd, sizeof(datum_dim));

  /* check that the dimension parameter array was allocated. */
  if (D->dims == NULL)
    throw("failed to allocate %u datum dimensions", D->nd);

  /* store the dimension information. */
  for (d = 0; d < D->nd; d++) {
    /* store the size parameters. */
    D->dims[d].sz = D->dims[d].td = D->dims[d].tdunif =
      (unsigned int) hdr.dims[d].sz;

    /* store the status flags. */
    D->dims[d].ft = 1;

    /* store the spectral parameters. */
    D->dims[d].carrier = (real) hdr.dims[d].sf;
    D->dims[d].width = (real) hdr.dims[d].sw;
    D->dims[d].offset = (real) hdr.dims[d].ref;

    /* check if the offset value needs to be converted. */
    switch (hdr.dims[d].refunits) {
      /* ppm -> hz. */
      case NV_REFUNIT_PPM:
        D->dims[d].offset *= D->dims[d].carrier;
        break;

      /* other. */
      default:
        break;
    }

    /* store the nucleus string. */
    strncpy(D->dims[d].nuc, hdr.dims[d].label, 8);
    D->dims[d].nuc[7] = '\0';
  }

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
  /* FIXME: implement nv_fwrite_datum() */throw("unimplemented!");

  /* return success. */
  return 1;
}

