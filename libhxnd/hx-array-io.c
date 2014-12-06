
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

/* include the n-dimensional math header. */
#include <hxnd/hx.h>

/* include the byte-level data header. */
#include <hxnd/bytes.h>

/* define the number of (u64) members in the header of binary array files.
 */
#define HX_ARRAY_FWRITE_SZ_HDR  6

/* define the size of the data buffer for raw file reading, in bytes.
 */
#define HX_ARRAY_FREAD_SZ_BUF  33554432

/* hx_array_print(): prints a hypercomplex multidimensional array as text.
 * @x: the array to print data from.
 * @fname: the output filename.
 */
int hx_array_print (hx_array *x, const char *fname) {
  /* declare a few required variables:
   * @fh: the file handle used for writing.
   */
  int *arr, idx, i;
  FILE *fh;

  /* open the output file. */
  if (fname)
    fh = fopen(fname, "wb");
  else
    fh = stdout;

  /* check that the file was opened. */
  if (!fh)
    throw("failed to open '%s'", fname);

  /* allocate an index array for iteration. */
  arr = hx_array_index_alloc(x->k);

  /* check that allocation was successful. */
  if (!arr)
    throw("failed to allocate %d indices", x->k);

  /* iterate over the points of the array. */
  idx = 0;
  do {
    /* print the indices. */
    for (i = 0; i < x->k; i++)
      fprintf(fh, "%6d ", arr[i]);

    /* print the coefficients. */
    for (i = 0; i < x->n; i++)
      fprintf(fh, "%18.8e ", x->x[i + x->n * idx]);

    /* print a newline. */
    fprintf(fh, "\n");

    /* increment the linear index. */
    idx++;
  } while (hx_array_index_incr(x->k, x->sz, arr));

  /* free the index array. */
  free(arr);

  /* close the output file. */
  if (fname)
    fclose(fh);

  /* return success. */
  return 1;
}

/* hx_array_check_magic(): checks the first "magic" bytes of a file and
 * returns whether they match the hypercomplex array magic number.
 * @fname: the input filename.
 */
int hx_array_check_magic (const char *fname) {
  /* declare a few required variables:
   * @fh: input file handle.
   * @wd: first word of input file.
   */
  uint64_t wd;
  FILE *fh;

  /* open the input file. */
  fh = fopen(fname, "rb");

  /* check that the file was opened. */
  if (!fh)
    throw("failed to open '%s'", fname);

  /* read the first word. */
  if (!fread(&wd, sizeof(uint64_t), 1, fh))
    throw("failed to read magic number");

  /* close the input file. */
  fclose(fh);

  /* check the magic word, without swapping. */
  if (wd == HX_ARRAY_MAGIC)
    return 1;

  /* swap the word and check again. */
  bytes_swap_u64(&wd);
  if (wd == HX_ARRAY_MAGIC)
    return 1;

  /* no match. */
  return 0;
}

/* hx_array_fwrite(): writes a hypercomplex multidimensional array to an
 * opened file stream.
 * @x: pointer to the source array.
 * @fh: the output file stream.
 */
int hx_array_fwrite (hx_array *x, FILE *fh) {
  /* declare a few required variables:
   * @wd: header that contains all array properties.
   * @n_wd: number of words in the header.
   * @i: word index in the header.
   * @k: dimension index.
   */
  unsigned int n_wd;
  uint64_t *wd;
  int i, k;

  /* allocate memory for the header bytes. */
  n_wd = HX_ARRAY_FWRITE_SZ_HDR + x->k;
  wd = (uint64_t*) calloc(n_wd, sizeof(uint64_t));

  /* check that the array was allocated. */
  if (!wd)
    throw("failed to allocate %d-word header", n_wd);

  /* initialize the word index. */
  i = 0;

  /* store a magic number into the header: 'HXNDARRY' */
  wd[i++] = (uint64_t) HX_ARRAY_MAGIC;

  /* store the array properties into the header. */
  wd[i++] = (uint64_t) x->d;
  wd[i++] = (uint64_t) x->n;
  wd[i++] = (uint64_t) x->k;
  wd[i++] = (uint64_t) x->len;
  wd[i++] = (uint64_t) sizeof(real);

  /* store the array dimension sizes into the header. */
  for (k = 0; k < x->k; k++)
    wd[i++] = (uint64_t) x->sz[k];

  /* write the header. */
  if (fwrite(wd, sizeof(uint64_t), n_wd, fh) != n_wd)
    throw("failed to write %d header words", n_wd);

  /* write the array data. */
  if (fwrite(x->x, sizeof(real), x->len, fh) != x->len)
    throw("failed to write %d reals", x->len);

  /* free the header array. */
  free(wd);

  /* return success. */
  return 1;
}

/* hx_array_fread(): reads a hypercomplex multidimensional array from an
 * opened file stream.
 * @x: pointer to the destination array.
 * @fh: the input file stream.
 */
int hx_array_fread (hx_array *x, FILE *fh) {
  /* declare a few required variables:
   * @wd: header that contains all array properties.
   * @i: word index in the header.
   * @k: dimension index.
   */
  unsigned int swapping;
  uint64_t *wd1, wd0[HX_ARRAY_FWRITE_SZ_HDR];
  int i, k, n_read;

  /* read the first six words from the file. */
  n_read = fread(wd0, sizeof(uint64_t), HX_ARRAY_FWRITE_SZ_HDR, fh);
  if (n_read != HX_ARRAY_FWRITE_SZ_HDR)
    throw("failed to read initial header words");

  /* check the first word in the header. if it does not match, then
   * byte-swap the data.
   */
  if (wd0[0] != HX_ARRAY_MAGIC) {
    /* no match. swap the bytes of each word. */
    bytes_swap((uint8_t*) wd0, HX_ARRAY_FWRITE_SZ_HDR, sizeof(uint64_t));

    /* now check the magic word. */
    if (wd0[0] != HX_ARRAY_MAGIC)
      throw("invalid magic number 0x%016lx", wd0[0]);

    /* set the swapping flag. */
    swapping = 1;
  }
  else {
    /* match. no swaps needed. */
    swapping = 0;
  }

  /* initialize the word index. */
  i = 1;

  /* read the array parameters from the header. */
  x->d = (int) wd0[i++];
  x->n = (int) wd0[i++];
  x->k = (int) wd0[i++];
  x->len = (int) wd0[i++];

  /* check that the data type size matches ours. if not, fail miserably.
   */
  if (wd0[i++] != sizeof(real))
    throw("word size mismatch (%u != %u)", wd0[i - 1], sizeof(real));

  /* allocate memory for the rest of the header. */
  wd1 = (uint64_t*) calloc(x->k, sizeof(uint64_t));

  /* check that the array memory was allocated. */
  if (!wd1)
    throw("failed to allocate %d-word header", x->k);

  /* read the entire header from the file. */
  if (fread(wd1, sizeof(uint64_t), x->k, fh) != x->k)
    throw("failed to read %d header words", x->k);

  /* byte-swap, if required. */
  if (swapping)
    bytes_swap((uint8_t*) wd1, x->k, sizeof(uint64_t));

  /* allocate memory for the dimension sizes array. */
  x->sz = (int*) calloc(x->k, sizeof(int));

  /* ensure the allocation was successful. */
  if (x->sz == NULL)
    throw("failed to allocate %d sizes", x->k);

  /* read the dimension sizes from the header. */
  for (k = 0, i = 0; k < x->k; k++)
    x->sz[k] = wd1[i++];

  /* allocate memory for the array data. */
  x->x = (real*) calloc(x->len, sizeof(real));

  /* ensure the allocation was successful. */
  if (!x->x)
    throw("failed to allocate %d reals", x->len);

  /* read the array data from the file. */
  if (fread(x->x, sizeof(real), x->len, fh) != x->len)
    throw("failed to read %d reals", x->len);

  /* byte-swap, if required. */
  if (swapping)
    bytes_swap((uint8_t*) x->x, x->len, sizeof(real));

  /* ensure that the d-dimensional shared multiplication table has been
   * initialized, and return failure if not.
   */
  if (!(x->tbl = hx_algebras_get(x->d)))
    throw("failed to retrieve %d-algebra", x->d);

  /* free the header array. */
  free(wd1);

  /* and return success. */
  return 1;
}

/* hx_array_save(): saves a hypercomplex multidimensional array to a file,
 * or standard output if a NULL filename was passed.
 * @x: the array to save data from.
 * @fname: the output filename.
 */
int hx_array_save (hx_array *x, const char *fname) {
  /* declare a required variable:
   * @fh: file handle used for writing.
   */
  FILE *fh;

  /* open the output file. */
  if (fname)
    fh = fopen(fname, "wb");
  else
    fh = stdout;

  /* check that the file was opened. */
  if (!fh)
    throw("failed to open '%s'", fname);

  /* write the data to the file. */
  if (!hx_array_fwrite(x, fh))
    throw("failed to write '%s'", fname);

  /* close the output file. */
  if (fname)
    fclose(fh);

  /* return success. */
  return 1;
}

/* hx_array_load(): loads a hypercomplex multidimensional array from a file,
 * or standard input if a NULL filename is passed.
 * @x: the array to load data into.
 * @fname: the input filename.
 */
int hx_array_load (hx_array *x, const char *fname) {
  /* declare a required variable:
   * @fh: file handle used for writing.
   */
  FILE *fh;

  /* open the input file. */
  if (fname)
    fh = fopen(fname, "rb");
  else
    fh = stdin;

  /* check that the file was opened. */
  if (!fh)
    throw("failed to open '%s'", fname);

  /* read the contents of the file. */
  if (!hx_array_fread(x, fh))
    throw("failed to read '%s'", fname);

  /* close the input file. */
  if (fname)
    fclose(fh);

  /* return success. */
  return 1;
}

/* hx_array_fread_raw(): read raw bytes from a file out of a given format
 * and into a linear hypercomplex array.
 * @fh: input file handle.
 * @x: output array structure pointer.
 * @endian: raw word byte ordering.
 * @wordsz: number of bytes per word.
 * @isflt: whether the words are floats.
 * @offhead: file header offset in bytes.
 * @offblk: file block offset in bytes.
 * @nblks: number of data blocks.
 * @nwords: number of data words per block.
 * @nalign: number of bytes to align blocks to.
 */
int hx_array_fread_raw (FILE *fh, hx_array *x, enum byteorder endian,
                        unsigned int wordsz, unsigned int isflt,
                        unsigned int offhead, unsigned int offblk,
                        unsigned int nblks, unsigned int nwords,
                        unsigned int nalign) {
  /* declare variables required for buffered reading:
   * @pos: current file position, in bytes.
   * @i: main block loop counter.
   * @n: number of blocks to read.
   * @nbytes: number of bytes per block.
   * @nbuf: number of buffer bytes.
   * @nask: number to requested bytes.
   * @nread: number of read bytes.
   * @nrem: number of bytes remaining.
   * @ialign: number of alignment offsets.
   * @buf: array of buffer bytes.
   */
  unsigned int pos, i, n, nbytes, nbuf, nask, nread, nrem, ialign;
  uint8_t *buf;

  /* declare variables required for coefficient allocation and access:
   * @k: raw buffer byte index for conversion into array words.
   * @xi: current array coefficient pointer.
   * @len: total number of words.
   */
  unsigned int k;
  real *xi;
  int len;

  /* compute the number of bytes per block, and the total byte count. */
  nbytes = nwords * wordsz;
  nrem = nbytes * nblks;

  /* compute the total number of data words. */
  len = nblks * nwords;

  /* initialize the block count and size. */
  nbuf = nbytes;
  n = nblks;

  /* determine whether block alignment is specified. */
  if (nalign) {
    /* check for the special case in which the reads may be lumped:
     * 1. the blocks fall on the alignment boundaries.
     * 2. no block header exists.
     * 3. no file header exists.
     */
    if (offhead == 0 && offblk == 0 && nbytes % nalign == 0) {
      /* in this special case, we can increase the buffer size. */
      nbuf *= n;
      n = 1;
    }
  }

  /* determine whether logical blocking must be done. */
  if (n == 1 && nbuf > HX_ARRAY_FREAD_SZ_BUF) {
    /* logical blocking is needed. determine the number of logical blocks:
     *  - if the buffer is evenly divided by the max size, use that ratio.
     *  - if not, add one to the ratio, where the last read will be partial.
     */
    n = (nbuf % HX_ARRAY_FREAD_SZ_BUF ?
         nbuf / HX_ARRAY_FREAD_SZ_BUF + 1 :
         nbuf / HX_ARRAY_FREAD_SZ_BUF);

    /* use the maximum allowed read buffer size. */
    nbuf = HX_ARRAY_FREAD_SZ_BUF;
  }

  /* allocate memory for the buffer bytes. */
  buf = (uint8_t*) calloc(nbuf, sizeof(uint8_t));

  /* check that memory was allocated. */
  if (!buf)
    throw("failed to allocate %u buffer bytes", nbytes);

  /* allocate the destination array structure. */
  if (!hx_array_alloc(x, 0, 1, &len))
    throw("failed to allocate %d array coefficients", len);

  /* move past the file header. */
  if (offhead && fseek(fh, offhead, SEEK_SET))
    throw("failed to seek %u bytes past file header", offhead);

  /* initialize the file position and the alignment offset. */
  pos = offhead;
  ialign = 0;

  /* loop over the data blocks. */
  for (i = 0, xi = x->x; i < n; i++) {
    /* check if alignment is required. */
    if (nalign) {
      /* align to the next boundary. */
      ialign = pos / nalign;
      while (pos % nalign)
        pos = ++ialign * nalign;

      /* move to the new aligned location. */
      if (fseek(fh, pos, SEEK_SET))
        throw("failed to seek to %u-byte alignment boundary", nalign);
    }

    /* move past the block header. */
    pos += offblk;
    if (offblk && fseek(fh, offblk, SEEK_CUR))
      throw("failed to seek %u bytes past block header", offblk);

    /* determine how many bytes to ask for. */
    nask = (nrem > nbuf ? nbuf : nrem);

    /* read the current data block. */
    nread = fread(buf, sizeof(uint8_t), nask, fh);

    /* check that the read succeeded. */
    if (nread != nask)
      throw("failed to read data block #%u", i);

    /* subtract the read bytes from the bytes remaining. */
    nrem -= nread;

    /* check if byte swaps are required. */
    if (!bytes_native(endian) && wordsz > 1)
      bytes_swap(buf, nbuf / wordsz, wordsz);

    /* copy the read words into the final array. */
    for (k = 0; k < nread; k += wordsz)
      *(xi++) = bytes_toword(buf + k, wordsz, isflt);
  }

  /* free the allocated buffer memory. */
  free(buf);

  /* return success. */
  return 1;
}

