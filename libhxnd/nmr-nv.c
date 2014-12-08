
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
                    enum byteorder *endian,
                    struct nv_header *hdr) {
  /* declare a few required variables:
   */
  unsigned int i, n_bytes, n_words;
  uint8_t *bytes;

  /* read in the file header bytes. */
  n_bytes = sizeof(struct nv_header);
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
    *endian = bytes_get_nonnative();
  }
  else {
    /* same endianness. */
    *endian = bytes_get_native();
  }

  /* free the read bytes. */
  free(bytes);

  /* return success. */
  return 1;
}

/* nv_copy_header(): copy the contents of an nmrview file header from one
 * structure to another.
 * @dst: destination header structure pointer.
 * @src: source header structure pointer.
 */
void nv_copy_header (struct nv_header *dst, struct nv_header *src) {
  /* perform the raw data copy. */
  memcpy(dst, src, sizeof(struct nv_header));
}

/* nv_adjust_header(): adjust the dimension sizes of an nmrview file header
 * in an attempt to evenly them divide into the block sizes and to match
 * the total file size with the expected file size.
 * @fname: the input filename.
 * @hdr: the header to adjust.
 */
int nv_adjust_header (const char *fname, struct nv_header *hdr) {
  /* declare a few required variables:
   * @i, @j: general-purpose loop counters.
   * @n_actual: number of data words in the input file.
   * @n_calc: estimated value of @n_actual based on header values.
   */
  unsigned int i, j, n_actual, n_calc;

  /* get the file size. */
  n_actual = bytes_size(fname);

  /* check that the file size is nonzero. */
  if (!n_actual || n_actual <= sizeof(struct nv_header))
    throw("invalid file size of %u-bytes", n_actual);

  /* subtract the number of bytes in the header. */
  n_actual -= sizeof(struct nv_header);

  /* check that the remaining size is divided by the word size. */
  if (n_actual % sizeof(float))
    throw("invalid data size of %u bytes", n_actual);

  /* divide the data section size by the word size. */
  n_actual /= sizeof(float);

  /* compute the estimated number of data words. */
  for (i = 0, n_calc = 1; i < hdr->ndims; i++)
    n_calc *= hdr->dims[i].sz;

  /* check if the estimated size matches the true size. */
  if (n_calc == n_actual)
    return 1;

  /* adjust the header values of each dimension. */
  for (i = 0; i < hdr->ndims; i++) {
    /* skip dimensions that are already evenly block-divided. */
    if (hdr->dims[i].sz % hdr->dims[i].szblk == 0)
      continue;

    /* compute an adjusted number of points in the current dimension. */
    hdr->dims[i].sz = n_actual;
    for (j = 1; j < hdr->ndims; j++)
      hdr->dims[i].sz /= hdr->dims[(i + j) % hdr->ndims].sz;

    /* fail if the adjusted size does not even divide into blocks. */
    if (hdr->dims[i].sz % hdr->dims[i].szblk)
      throw("adjusted size %d (#%u) does not evenly divide block size %d",
            hdr->dims[i].sz, i, hdr->dims[i].szblk);
  }

  /* compute the estimated number of data words. */
  for (i = 0, n_calc = 1; i < hdr->ndims; i++)
    n_calc *= hdr->dims[i].sz;

  /* check if the estimated size matches the true size. */
  if (n_calc != n_actual)
    throw("adjustment failed to correct file size mismatch");

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
int nv_tiler (hx_array *x, struct nv_header *hdr, int dir) {
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
  if (!hx_array_tiler(x, k, nt, szt, dir, HX_ARRAY_INCR_NORMAL))
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
  struct nv_header hdr;

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

  /* adjust the file header, if necessary. */
  if (!nv_adjust_header(fname, &hdr))
    throw("failed to perform header adjustment");

  /* compute the byte offset from which to begin reading point data. */
  offset = sizeof(struct nv_header);

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

/* nv_fix_adjustment(): undo any size adjustments made to the array of a
 * datum structure in order to successfully load an nmrview-format file.
 * @D: pointer to the datum structure.
 */
int nv_fix_adjustment (datum *D) {
  /* declare a few required variables:
   * @d: dimension loop counter.
   * @adjusted: whether the sizes have been adjusted.
   * @sznew: array of new datum array sizes.
   */
  unsigned int d, adjusted;
  int *sznew;

  /* determine whether adjustment was performed. */
  for (d = 0, adjusted = 0; d < D->nd; d++) {
    /* if any single dimension contains differing @sz and @td values,
     * then an adjustment was performed.
     */
    if (D->dims[d].sz != D->dims[d].td)
      adjusted = 1;
  }

  /* return successfully if no adjustments were made. */
  if (!adjusted)
    return 1;

  /* allocate the new size array. */
  sznew = hx_array_index_alloc(D->array.k);

  /* check that the size array was allocated. */
  if (!sznew)
    throw("failed to allocate %d indices", D->array.k);

  /* fill the new size array. */
  for (d = 0; d < D->nd; d++) {
    /* restore the size value and store its value. */
    D->dims[d].sz = D->dims[d].td;
    sznew[d] = (int) D->dims[d].sz;
  }

  /* resize the core array. */
  if (!hx_array_resize(&D->array, D->array.d, D->array.k, sznew))
    throw("failed to resize core datum array");

  /* free the allocated size array. */
  free(sznew);

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
   * @endian: the byte ordering of the data file.
   * @hdr: the original file header structure.
   * @hdradj: the adjusted file header structure.
   * @d: dimension loop counter.
   * @ts: calendar time structure to hold header values.
   */
  enum byteorder endian = BYTES_ENDIAN_AUTO;
  struct nv_header hdr, hdradj;
  unsigned int d;
  struct tm ts;

  /* read the header information from the data file. */
  if (!nv_read_header(fname, &endian, &hdr))
    throw("failed to read header of '%s'", fname);

  /* copy the header information into the structure to be adjusted. */
  nv_copy_header(&hdradj, &hdr);

  /* adjust the file header, if necessary. */
  if (!nv_adjust_header(fname, &hdradj))
    throw("failed to perform header adjustment");

  /* check if time/date information exists in the header fields. */
  if (hdr.year && hdr.month) {
    /* parse the date header fields. */
    ts.tm_mon = (int) hdr.month - 1;
    ts.tm_mday = (int) hdr.day;
    ts.tm_year = (int) hdr.year - 1900;

    /* convert the date and time into an epoch offset. */
    D->epoch = mktime(&ts);
  }

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
    /* store the adjusted size parameter. this will allow for successful
     * refactoring of the array content.
     */
    D->dims[d].sz = (unsigned int) hdradj.dims[d].sz;

    /* store the unadjusted size parameter. this will be used by
     * nv_deadjust_datum() to crop the array of files that required
     * header adjustment.
     */
    D->dims[d].td = D->dims[d].tdunif = (unsigned int) hdr.dims[d].sz;

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
  D->endian = endian;

  /* return success. */
  return 1;
}

/* nv_fwrite_datum(): write an NMR datum structure in nmrview-format to
 * an opened file stream.
 * @D: pointer to the source structure.
 * @fh: the output file stream.
 */
int nv_fwrite_datum (datum *D, FILE *fh) {
  /* declare variables required to write header information:
   * @d: dimension loop counter.
   * @i_div: current division to subdivide.
   * @n_tile: number of words per tile.
   * @hdr: the output file header structure.
   * @ts: calendar time structure.
   */
  unsigned int i_div, n_tile;
  struct nv_header hdr;
  struct tm *ts;

  /* declare variables required for size adjustment:
   * @szadj: whether to adjust the array size.
   * @sz: previous array size.
   * @sznew: adjusted array size.
   */
  int szadj, *sz, *sznew;

  /* declare a dimension loop counter. */
  unsigned int d;

  /* declare an array structure to use for writing. */
  hx_array xout;

  /* allocate size arrays for possible adjustments. */
  sz = hx_array_index_alloc(D->array.k);
  sznew = hx_array_index_alloc(D->array.k);

  /* check that array allocation succeeded. */
  if (!sz || !sznew)
    throw("failed to allocate two sets of %d indices", D->array.k);

  /* initialize the file header. */
  memset(&hdr, 0, sizeof(struct nv_header));

  /* set the fields of the file header. */
  hdr.magic = NV_MAGIC;
  hdr.fhdrsz = sizeof(struct nv_header);
  hdr.bhdrsz = 0;
  hdr.ndims = D->nd;

  /* compute the calendar date structure fields. */
  ts = gmtime(&D->epoch);

  /* store the header date fields. */
  hdr.month = ts->tm_mon + 1;
  hdr.day = ts->tm_mday;
  hdr.year = ts->tm_year + 1900;

  /* set the extraneous file header fields. */
  hdr.temp = 298.15;
  strcpy(hdr.sequence, "");
  strcpy(hdr.comment, "");

  /* initialize the tile sizes. */
  for (d = 0; d < D->nd; d++)
    hdr.dims[d].szblk = 2;

  /* determine the tile size of each dimension. */
  i_div = 0;
  do {
    /* multiply the tile size of the currently indexed dimension. */
    hdr.dims[i_div].szblk *= 2;

    /* increment the dimension index. */
    if (++i_div >= D->nd)
      i_div = 0;

    /* compute the new tile size. */
    for (d = 0, n_tile = 2; d < D->nd; d++)
      n_tile *= hdr.dims[d].szblk;
  } while (n_tile < NV_MAX_TILE);

  /* compute the number of block elements. */
  for (d = 0, hdr.blkelem = 1; d < D->nd; d++)
    hdr.blkelem *= hdr.dims[d].szblk;

  /* configure each dimension sub-header. */
  for (d = 0, szadj = 0; d < D->nd; d++) {
    /* store the old dimension size. */
    sznew[d] = sz[d] = D->dims[d].sz;

    /* if any single dimension does not evenly divide into blocks,
     * the array will need to be resized.
     */
    if (sznew[d] % hdr.dims[d].szblk) {
      /* adjustment is required. */
      szadj = 1;

      /* compute the adjusted dimension size. */
      sznew[d] = hdr.dims[d].szblk * (sz[d] / hdr.dims[d].szblk + 1);
    }

    /* set the dimension sub-header fields. */
    hdr.dims[d].sz = sznew[d];
    hdr.dims[d].nblk = 16;
    hdr.dims[d].offblk = (d > 0 ? d * 16 : 1);
    hdr.dims[d].maskblk = hdr.dims[d].szblk - 1;
    hdr.dims[d].ptoff = 1 + d * 4;

    /* set the dimension spectral parameters. */
    hdr.dims[d].sf = D->dims[d].carrier;
    hdr.dims[d].sw = D->dims[d].width;
    hdr.dims[d].refpt = hdr.dims[d].sz / 2;
    hdr.dims[d].ref = D->dims[d].offset / D->dims[d].carrier;
    hdr.dims[d].refunits = NV_REFUNIT_PPM;

    /* set the spectral folding parameters. */
    hdr.dims[d].foldup = 0.0;
    hdr.dims[d].folddown = 0.0;

    /* copy the dimension label string. */
    strcpy(hdr.dims[d].label, D->dims[d].nuc);
  }

  /* check if the datum array is real. */
  if (D->array.d == 0) {
    /* just use the datum array. */
    if (!hx_array_copy(&xout, &D->array))
      throw("failed to copy core array");
  }
  else {
    /* copy the real component of the datum array. */
    if (!hx_array_copy_real(&xout, &D->array))
      throw("failed to copy real component of core array");
  }

  /* check if size adjustment is required. */
  if (szadj && !hx_array_resize(&xout, xout.d, xout.k, sznew))
    throw("failed to adjust array size for output");

  /* map the linear datum array back into tiles. */
  if (!nv_tileize(&xout, &hdr))
    throw("failed to delinearize array into tiles");

  /* de-adjust the header size values. */
  for (d = 0; d < D->nd; d++)
    hdr.dims[d].sz = sz[d];

  /* write the file header. */
  if (fwrite(&hdr, sizeof(struct nv_header), 1, fh) != 1)
    throw("failed to write file header");

  /* write the array data. */
  if (!hx_array_fwrite_raw(fh, &xout, bytes_get_native(), sizeof(float), 1))
    throw("failed to write core array data");

  /* throw away the real copy of the datum array. */
  hx_array_free(&xout);

  /* free the allocated size arrays. */
  free(sz);
  free(sznew);

  /* return success. */
  return 1;
}

