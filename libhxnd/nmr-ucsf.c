
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

/* include the ucsf header. */
#include <hxnd/nmr-ucsf.h>

/* define constants which determine the behavior of ucsf_map_tiles().
 */
#define UCSF_TILEMAP_FORWARD  (+1)
#define UCSF_TILEMAP_REVERSE  (-1)

/* ucsf_check_magic(): checks the first "magic" bytes of a file and
 * returns whether they match the ucsf-format file magic number.
 * @fname: the input filename.
 */
int ucsf_check_magic (const char *fname) {
  /* declare a few required variables:
   * @fh: input file handle.
   * @wd: first words of input file.
   */
  char wd[UCSF_NUM_MAGIC];
  FILE *fh;

  /* open the input file. */
  fh = fopen(fname, "rb");

  /* check that the file was opened. */
  if (!fh)
    throw("failed to open '%s'", fname);

  /* read the first ten words. */
  if (!fread(&wd, sizeof(char), UCSF_NUM_MAGIC, fh))
    throw("failed to read magic numbers");

  /* close the input file. */
  fclose(fh);

  /* check the magic string, without swapping. */
  if (strncmp(wd, UCSF_MAGIC, UCSF_NUM_MAGIC) == 0)
    return 1;

  /* no match. */
  return 0;
}

/* ucsf_read_header(): read the contents of a header from a ucsf-format file.
 * @fname: the input filename.
 * @endianness: byte order result pointer.
 * @hdr: file header result pointer.
 * @dims: dimension header array result pointer.
 */
int ucsf_read_header (const char *fname, enum byteorder *endianness,
                      struct ucsf_file_header *hdr,
                      struct ucsf_dim_header **dims) {
  /* declare a few required variables:
   */
  struct ucsf_dim_header *ldims;
  unsigned int n_total, n_calc, n_hdr, n_dims, i;
  uint8_t *bytes;

  /* read in the total number of bytes in the file. */
  n_total = bytes_size(fname);

  if (!n_total)
    throw("failed to read size of '%s'", fname);

  /* read in the file header bytes. */
  n_hdr = sizeof(struct ucsf_file_header);
  bytes = bytes_read_block(fname, 0, n_hdr);

  /* check that the bytes were read successfully. */
  if (!bytes)
    throw("failed to read file header from '%s'", fname);

  /* copy the header bytes onto the header structure. */
  memcpy(hdr, bytes, n_hdr);

  /* free the read bytes. */
  free(bytes);

  /* initialize the number of bytes to read per dimension header. */
  n_dims = sizeof(struct ucsf_dim_header);

  /* allocate the array of dimension headers. */
  ldims = *dims = (struct ucsf_dim_header*) malloc(hdr->ndims * n_dims);

  /* check that allocation succeeded. */
  if (ldims == NULL)
    throw("failed to allocate %u dimension headers", hdr->ndims);

  /* read in the dimension header bytes. */
  for (i = 0; i < hdr->ndims; i++) {
    /* read in the dimension header bytes. */
    bytes = bytes_read_block(fname, n_hdr + i * n_dims, n_dims);

    /* check that the bytes were read successfully. */
    if (!bytes)
      throw("failed to read dimension header %u from '%s'", i, fname);

    /* copy the header bytes onto the header structure. */
    memcpy(&ldims[i], bytes, n_dims);

    /* free the read bytes. */
    free(bytes);
  }

  /* calculate the number of expected total points. */
  for (i = 0, n_calc = sizeof(float); i < hdr->ndims; i++)
    n_calc *= ldims[i].npts;

  /* calculate the number of expected total bytes. */
  n_calc += n_hdr + hdr->ndims * n_dims;

  /* check if the calculated size matches the true size. */
  if (n_calc != n_total) {
    /* loop through the dimensions and swap the bytes. */
    for (i = 0, n_calc = sizeof(float); i < hdr->ndims; i++) {
      /* swap bytes of the header entries. */
      bytes_swap_u32(&ldims[i].npts);
      bytes_swap_u32(&ldims[i].sztile);
      bytes_swap_u32((uint32_t*) &ldims[i].carrier);
      bytes_swap_u32((uint32_t*) &ldims[i].width);
      bytes_swap_u32((uint32_t*) &ldims[i].center);

      /* compute the new calculated size. */
      n_calc *= ldims[i].npts;
    }

    /* calculate the number of expected total bytes. */
    n_calc += n_hdr + hdr->ndims * n_dims;

    /* check that the calculated size matches after swapping. */
    if (n_calc != n_total)
      throw("invalid file size of %u bytes", n_total);

    /* opposite endianness. */
    *endianness = bytes_get_nonnative();
  }
  else {
    /* same endianness. */
    *endianness = bytes_get_native();
  }

  /* return success. */
  return 1;
}

/* ucsf_map_tiles(): redoes or undoes the tiling inherent in ucsf-format
 * files, effectively mapping between tiled array data and a real linear
 * array suitable for datum_refactor_array().
 * @x: pointer to the array to linearize.
 * @fhdr: pointer to the file header.
 * @dhdr: array of dimension headers.
 * @dir: direction, either +1 (linearize) or -1 (tileize).
 */
int ucsf_map_tiles (hx_array *x,
                    struct ucsf_file_header *fhdr,
                    struct ucsf_dim_header *dhdr,
                    int dir) {
  /* declare a few required variables:
   * @i: general purpose loop counter.
   * @k: number of dimensions.
   * @idxi: input array linear index.
   * @idxo: output array linear index.
   * @sz: array of tile sizes.
   * @szt: array of tile counts.
   * @arr: multidimensional point index array.
   * @arrt: multidimensional tile index array.
   * @xcpy: temporary tile array.
   */
  int i, k, n, ncpy, idxi, idxo, *sz, *szt, *arr, *arrt;
  hx_array xcpy;

  /* gain a handle on the dimensionality of the array data. */
  k = (int) fhdr->ndims;

  /* allocate the size index arrays. */
  sz = hx_array_index_alloc(k);
  szt = hx_array_index_alloc(k);

  /* allocate the loop index arrays. */
  arr = hx_array_index_alloc(k);
  arrt = hx_array_index_alloc(k);

  /* check that all index allocations were successful. */
  if (!sz || !szt || !arr || !arrt)
    throw("failed to allocate four sets of %d indices", k);

  /* initialize the size arrays. */
  for (i = 0; i < k; i++) {
    /* check that the point count is evenly divided by the tile size. */
    if (dhdr[i].npts % dhdr[i].sztile)
      throw("tile size %u does not evenly divide point count %u",
            dhdr[i].sztile, dhdr[i].npts);

    /* set the tile point count. */
    sz[i] = dhdr[i].sztile;

    /* set the tile count. */
    szt[i] = dhdr[i].npts / dhdr[i].sztile;
  }

  /* allocate a temporary array for storing tile data. */
  if (!hx_array_copy(&xcpy, x))
    throw("failed to allocate tile array");

  /* compute the coefficient count and byte count per scalar.
   */
  n = x->n;
  ncpy = n * sizeof(real);

  /* loop over the data tiles. initialize the input linear index. */
  idxi = 0;
  do {
    /* loop over the tile points. */
    do {
      /* pack the tiled indices into a linear index. */
      hx_array_index_pack_tiled(k, szt, sz, arr, arrt, &idxo);

      /* determine the copy direction. */
      switch (dir) {
        /* forward: tiles -> linear. */
        case UCSF_TILEMAP_FORWARD:
          /* copy coefficients from the tiled array into the linear array. */
          memcpy(x->x + idxo * n, xcpy.x + idxi * n, ncpy);
          break;

        /* reverse: linear -> tiles. */
        case UCSF_TILEMAP_REVERSE:
          /* copy coefficients from the linear array into the tiled array. */
          memcpy(x->x + idxi * n, xcpy.x + idxo * n, ncpy);
          break;
      }

      /* increment the tile linear index. */
      idxi++;
    } while (hx_array_index_incr_rev(k, sz, arr));
  } while (hx_array_index_incr_rev(k, szt, arrt));

  /* free the allocated tile array. */
  hx_array_free(&xcpy);

  /* free the allocated index arrays. */
  free(sz);
  free(szt);
  free(arr);
  free(arrt);

  /* return success. */
  return 1;
}

/* ucsf_linearize(): maps tiles to linear array data.
 */
#define ucsf_linearize(x, fhdr, dhdr) \
  ucsf_map_tiles(x, fhdr, dhdr, UCSF_TILEMAP_FORWARD)

/* ucsf_delinearize(): maps linear array data to tiles.
 */
#define ucsf_delinearize(x, fhdr, dhdr) \
  ucsf_map_tiles(x, fhdr, dhdr, UCSF_TILEMAP_REVERSE)

/* ucsf_read(): reads a ucsf-format data file into a real linear array.
 * @fname: the input data filename.
 * @x: the output array.
 */
int ucsf_read (const char *fname, hx_array *x) {
  /* declare variables required to read headers.
   * @endianness: byte ordering of the input file.
   * @fhdr: file header structure.
   * @dhdr: array of dimension header structures.
   */
  enum byteorder endianness = BYTES_ENDIAN_AUTO;
  struct ucsf_file_header fhdr;
  struct ucsf_dim_header *dhdr;

  /* declare variables for reading raw data bytes:
   * @offset: byte offset where point data begins.
   * @i: general purpose loop counter.
   * @n: number of bytes to read.
   */
  unsigned int offset, i, n;
  uint8_t *bytes;

  /* attempt to read the file and dimension headers from the file. */
  if (!ucsf_read_header(fname, &endianness, &fhdr, &dhdr))
    throw("failed to read header of '%s'", fname);

  /* compute the byte offset from which to begin reading point data. */
  offset = sizeof(struct ucsf_file_header);
  offset += fhdr.ndims * sizeof(struct ucsf_dim_header);

  /* compute the number of bytes to read. */
  for (i = 0, n = sizeof(float); i < fhdr.ndims; i++)
    n *= dhdr[i].npts;

  /* read the data bytes in from the file. */
  bytes = bytes_read_block(fname, offset, n);

  /* check that the bytes were read successfully. */
  if (!bytes)
    throw("failed to read %u bytes from '%s'", n, fname);

  /* build a real linear array from the byte data. */
  if (!bytes_toarray(bytes, n, endianness, sizeof(float), 1, x))
    throw("failed to convert bytes to array");

  /* free the read byte data. */
  free(bytes);

  /* convert tiles to linear values. */
  if (!ucsf_linearize(x, &fhdr, dhdr))
    throw("failed to linearize tiled array");

  /* free the dimension header array. */
  free(dhdr);

  /* return success. */
  return 1;
}

/* ucsf_fill_datum(): intelligently parses ucsf parameters into an
 * NMR datum structure.
 * @fname: the input filename.
 * @D: pointer to the datum struct to fill.
 */
int ucsf_fill_datum (const char *fname, datum *D) {
  /* declare variables required to read headers:
   */
  enum byteorder endianness = BYTES_ENDIAN_AUTO;
  struct ucsf_file_header fhdr;
  struct ucsf_dim_header *dhdr;
  unsigned int d;

  /* attempt to read the file and dimension headers from the file. */
  if (!ucsf_read_header(fname, &endianness, &fhdr, &dhdr))
    throw("failed to read header of '%s'", fname);

  /* check the component count. */
  if (fhdr.ncomp != 1)
    throw("invalid data component count %u", fhdr.ncomp);

  /* store the dimensionality. */
  D->nd = (unsigned int) fhdr.ndims;

  /* check the dimensionality. */
  if ((int) fhdr.ndims < 1)
    throw("invalid dimensionality %u", fhdr.ndims);

  /* allocate the dimension parameter array. */
  D->dims = (datum_dim*) calloc(D->nd, sizeof(datum_dim));

  /* check that the dimension parameter array was allocated. */
  if (D->dims == NULL)
    throw("failed to allocate %u datum dimensions", D->nd);

  /* store the dimension information. */
  for (d = 0; d < D->nd; d++) {
    /* store the size parameters. */
    D->dims[d].sz = D->dims[d].td = D->dims[d].tdunif =
      (unsigned int) dhdr[d].npts;

    /* store the status flags. */
    D->dims[d].ft = 1;

    /* store the spectral parameters. */
    D->dims[d].carrier = (real) dhdr[d].carrier;
    D->dims[d].width = (real) dhdr[d].width;
    D->dims[d].offset = ((real) dhdr[d].center) * D->dims[d].carrier;

    /* store the nucleus string. */
    strncpy(D->dims[d].nuc, dhdr[d].nuc, 6);
    D->dims[d].nuc[7] = '\0';
  }

  /* free the dimension header array. */
  free(dhdr);

  /* store the filename string. */
  D->fname = (char*) malloc((strlen(fname) + 1) * sizeof(char));
  if (D->fname)
    strcpy(D->fname, fname);

  /* store the datum type. */
  D->type = DATUM_TYPE_UCSF;
  D->endian = endianness;

  /* return success. */
  return 1;
}

/* ucsf_fwrite_datum(): writes an NMR datum structure in ucsf-format to
 * an opened file stream.
 * @D: pointer to the source structure.
 * @fh: the output file stream.
 */
int ucsf_fwrite_datum (datum *D, FILE *fh) {
  /* declare variables required for header output.
   * @fhdr: output file header.
   * @dhdr: array of dimension headers.
   * @i_div: current dimension to subdivide.
   * @n_tile: number of bytes per tile.
   */
  struct ucsf_file_header fhdr;
  struct ucsf_dim_header *dhdr;
  unsigned int i_div, n_tile;

  /* declare a dimension loop counter. */
  unsigned int d;

  /* FIXME: handle writing of complex datum arrays. */

  /* allocate an array of dimension headers. */
  dhdr = (struct ucsf_dim_header*)
    calloc(D->nd, sizeof(struct ucsf_dim_header));

  /* check that allocation succeeded. */
  if (!dhdr)
    throw("failed to allocate %u dimension headers", D->nd);

  /* initialize the file header. */
  memset(&fhdr, 0, sizeof(struct ucsf_file_header));

  /* set the fields of the file header. */
  strcpy(fhdr.ftype, UCSF_MAGIC);
  fhdr.ndims = (uint8_t) D->nd;
  fhdr.ncomp = 1;
  fhdr.fmtver = 2;

  /* configure each dimension header. */
  for (d = 0; d < D->nd; d++) {
    /* set the nucleus string of the dimension header. */
    strcpy(dhdr[d].nuc, D->dims[d].nuc);

    /* set the point count of the dimension header. */
    dhdr[d].npts = (uint32_t) D->dims[d].sz;
    dhdr[d].sztile = dhdr[d].npts;

    /* set the spectral parameters of the dimension header. */
    dhdr[d].carrier = (float) D->dims[d].carrier;
    dhdr[d].width = (float) D->dims[d].width;
    dhdr[d].center = (float) D->dims[d].offset / dhdr[d].carrier;
  }

  /* check for zero-size dimensions. */
  for (d = 0; d < D->nd; d++) {
    /* fail if the dimension has zero size. */
    if (dhdr[d].npts == 0)
      throw("dimension %u has zero size", d);
  }

  /* determine the tile size of each dimension. */
  i_div = 0;
  do {
    /* do not attempt to subdivide odd tile sizes further. */
    while (dhdr[i_div].sztile % 2 && i_div < D->nd)
      i_div++;

    /* check that we found an index to subdivide. */
    if (i_div >= D->nd)
      throw("failed to identify suitable tile sizes");

    /* divide tile size of the currently indexed dimension. */
    dhdr[i_div].sztile /= 2;

    /* increment the dimension index. */
    if (++i_div >= D->nd)
      i_div = 0;

    /* compute the new tile size. */
    for (d = 0, n_tile = sizeof(float); d < D->nd; d++)
      n_tile *= dhdr[d].sztile;
  } while (n_tile > 32768);

  /* map the linear datum array back into tiles. */
  if (!ucsf_delinearize(&D->array, &fhdr, dhdr))
    throw("failed to delinearize array into tiles");

  /* write the file header. */
  if (fwrite(&fhdr, sizeof(struct ucsf_file_header), 1, fh) != 1)
    throw("failed to write file header");

  /* write the dimension headers. */
  if (fwrite(dhdr, sizeof(struct ucsf_dim_header), D->nd, fh) != D->nd)
    throw("failed to write %u dimension headers", D->nd);

  /* write the array data. */
  if (fwrite(D->array.x, sizeof(real), D->array.len, fh) != D->array.len)
    throw("failed to write core array data");

  /* map the tiled array back into the proper format, just in case it
   * needs to be used again.
   */
  if (!ucsf_linearize(&D->array, &fhdr, dhdr))
    throw("failed to re-linearize tiled array");

  /* free the dimension header array. */
  free(dhdr);

  /* return success. */
  return 1;
}

