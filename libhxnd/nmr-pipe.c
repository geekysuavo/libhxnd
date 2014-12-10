
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

/* define a magic number use to check byte ordering of pipe files. */
#define PIPE_MAGIC  2.345

/* define the maximum number of dimensions supported by pipe files. */
#define PIPE_MAXDIM  4

/* define sizes of the string struct members of the pipe header.
 */
#define PIPE_HDRSTR_SZ_LABEL      8
#define PIPE_HDRSTR_SZ_SRC       16
#define PIPE_HDRSTR_SZ_USER      16
#define PIPE_HDRSTR_SZ_OPER      32
#define PIPE_HDRSTR_SZ_COMMENT  160
#define PIPE_HDRSTR_SZ_TITLE     60

/* define the value of the @phase2d header for writing pipe files.
 */
#define PIPE_PHASE2D_STATES  2

/* define values that the @quad header variables may take in pipe files.
 */
#define PIPE_QUAD_COMPLEX  0
#define PIPE_QUAD_REAL     1
#define PIPE_QUAD_PSEUDO   2
#define PIPE_QUAD_SE       3
#define PIPE_QUAD_GE       4

/* define values that the @aqsgn header variables may take in pipe files.
 */
#define PIPE_AQSGN_NONE  0
#define PIPE_AQSGN_SEQ   1
#define PIPE_AQSGN_ALT   2
#define PIPE_AQSGN_NEG  16

/* pipe_header: properly padded header of data contained in a pipe-format file.
 * pipe file headers are packed as 512 32-bit floating point numbers, some of
 * which are, in fact, characters.
 */
struct pipe_header {
  /* (0) @magic: should be zero in valid nmrpipe data.
   * (1) @format: constant defining floating point format.
   * (2) @order: constant defining byte order in floats.
   */
  float magic;
  float format;
  float order;

  /* (3..8) @pad0 */
  float pad0[6];

  /* (9) @ndims: number of dimensions.
   * (10) @obs_f3: third-dimension observe frequency.
   * (11) @sw_f3: third-dimension spectral witdh.
   * (12) @orig_f3: third-dimension spectral origin.
   * (13) @ftflag_f3: third-dimension fft status.
   */
  float ndims;
  float obs_f3;
  float sw_f3;
  float orig_f3;
  float ftflag_f3;

  /* (14) @plane_loc: location of this plane.
   * (15) @size_f3: third-dimension size.
   */
  float plane_loc;
  float size_f3;

  /* (16..17) @label_f2: second-dimension label.
   * (18..19) @label_f1: first-dimension label.
   * (20..21) @label_f3: third-dimension label.
   * (22..23) @label_f4: fourth-dimension label.
   */
  char label_f2[PIPE_HDRSTR_SZ_LABEL];
  char label_f1[PIPE_HDRSTR_SZ_LABEL];
  char label_f3[PIPE_HDRSTR_SZ_LABEL];
  char label_f4[PIPE_HDRSTR_SZ_LABEL];

  /* (24..27) @dimorder: dimension ordering array. */
  float dimorder[4];

  /* (28) @obs_f4: fourth-dimension observe frequency.
   * (29) @sw_f4: fourth-dimension spectral witdh.
   * (30) @orig_f4: fourth-dimension spectral origin.
   * (31) @ftflag_f4: fourth-dimension fft status.
   * (32) @size_f4: fourth-dimension size.
   */
  float obs_f4;
  float sw_f4;
  float orig_f4;
  float ftflag_f4;
  float size_f4;

  /* (33..49) @pad1 */
  float pad1[17];

  /* (50) @apod_f3: third-dimension apodization point count.
   * (51) @quad_f3: third-dimension quadrature flag.
   */
  float apod_f3;
  float quad_f3;

  /* (52) @pad2 */
  float pad2;

  /* (53) @apod_f4: fourth-dimension apodization point count.
   * (54) @quad_f4: fourth-dimension quadrature flag.
   * (55) @quad_f1: first-dimension quadrature flag.
   * (56) @quad_f2: second-dimension quadrature flag.
   */
  float apod_f4;
  float quad_f4;
  float quad_f1;
  float quad_f2;

  /* (57) @pipe: whether the file is in a pipe stream.
   * (58) @units_f3: third-dimension axis unit type.
   * (59) @units_f4: fourth-dimension axis unit type.
   * (60) @ph0_f3: third-dimension zero-order phase correction.
   * (61) @ph1_f3: third-dimension first-order phase correction.
   * (62) @ph0_f4: fourth-dimension zero-order phase correction.
   * (63) @ph1_f4: fourth-dimension first-order phase correction.
   * (64) @aqsgn_f2: second-dimension fft sign adjustment.
   */
  float pipe;
  float units_f3;
  float units_f4;
  float ph0_f3;
  float ph1_f3;
  float ph0_f4;
  float ph1_f4;
  float aqsgn_f2;

  /* (65) @partition: slice count used in client-server mode. */
  float partition;

  /* (66) @car_f2: second-dimension carrier frequency.
   * (67) @car_f1: first-dimension carrier frequency.
   * (68) @car_f3: third-dimension carrier frequency.
   * (69) @car_f4: fourth-dimension carrier frequency.
   */
  float car_f2;
  float car_f1;
  float car_f3;
  float car_f4;

  /* (70..74) @user: general-purpose user data.
   * (75) @pipecount: number of functions in the pipe.
   * (76) @pad3
   */
  float user[5];
  float pipecount;
  float pad3;

  /* (77) @first_plane: first z-plane in the subset.
   * (78) @last_plane: last z-plane in the subset.
   */
  float first_plane;
  float last_plane;

  /* (79) @center_f2: second-dimension center frequency.
   * (80) @center_f1: first-dimension center frequency.
   * (81) @center_f3: third-dimension center frequency.
   * (82) @center_f4: fourth-dimension center frequency.
   */
  float center_f2;
  float center_f1;
  float center_f3;
  float center_f4;

  /* (83..94) @pad4 */
  float pad4[12];

  /* (95) @apod_f2: second-dimension apodization point count.
   * (96) @ftsz_f2: second-dimension fft points count.
   * (97) @realsz: number of valid time-domain points (obsolete).
   * (98) @ftsz_f1: first-dimension fft points count.
   * (99) @sz: number of points in the current dimension (real|imag).
   */
  float apod_f2;
  float ftsz_f2;
  float realsz;
  float ftsz_f1;
  float sz;

  /* (100) @sw_f2: second-dimension spectral witdh.
   * (101) @orig_f2: second-dimension spectral origin.
   */
  float sw_f2;
  float orig_f2;

  /* (102..105) @pad5 */
  float pad5[4];

  /* (106) @quad: quadrature flag.
   * (107) @pad6
   * (108) @zf_f2: second-dimension zero-fill value.
   * (109) @ph0_f2: second-dimension zero-order phase correction.
   * (110) @ph1_f2: second-dimension first-order phase correction.
   * (111) @lb_f2: second-dimension line-broadening value.
   */
  float quad;
  float pad6;
  float zf_f2;
  float ph0_f2;
  float ph1_f2;
  float lb_f2;

  /* (112..118) @pad7 */
  float pad7[7];

  /* (119) @obs_f2: second-dimension observe frequency. */
  float obs_f2;

  /* (120..134) @pad8 */
  float pad8[15];

  /* (135) @mcflag: magnitude calculation performed. */
  float mcflag;

  /* (136..151) @pad9 */
  float pad9[16];

  /* (152) @units_f2: second-dimension axis unit type.
   * (153) @noise: estimated noise floor intensity.
   */
  float units_f2;
  float noise;

  /* (154..156) @pad10 */
  float pad10[3];

  /* (157) @temperature: sample temperature in degrees celsius. */
  float temperature;

  /* (158..179) @pad11 */
  float pad11[22];

  /* (180) @rank: estimate of matrix rank. */
  float rank;

  /* (181..198) @pad12 */
  float pad12[18];

  /* (199) @tau: "a tau value for spectral series." (?)
   * (200) @ftsz_f3: third-dimension fft points count.
   * (201) @ftsz_f4: fourth-dimension fft points count.
   */
  float tau;
  float ftsz_f3;
  float ftsz_f4;

  /* (202..217) @pad13 */
  float pad13[16];

  /* (218) @obs_f1: first-dimension observe frequency.
   * (219) @specnum: number of complex 1d slices in the file.
   * (220) @ftflag_f2: second-dimension fft status.
   * (221) @trans: transposition status.
   * (222) @ftflag_f1: first-dimension fft status.
   */
  float obs_f1;
  float specnum;
  float ftflag_f2;
  float trans;
  float ftflag_f1;

  /* (223..228) @pad14 */
  float pad14[6];

  /* (229) @sw_f1: first-dimension spectral width. */
  float sw_f1;

  /* (230..233) @pad15 */
  float pad15[4];

  /* (234) @units_f1: first-dimension axis unit type. */
  float units_f1;

  /* (235..242) @pad16 */
  float pad16[8];

  /* (243) @lb_f1: first-dimension line-broadening value.
   * (244) @pad17
   * (245) @ph0_f1: first-dimension zero-order phase correction.
   * (246) @ph1_f1: second-dimension first-order phase correction.
   * (247) @max: maximum real component of the data.
   * (248) @min: minimum real component of the data.
   * (249) @orig_f1: first-dimension spectral origin.
   * (250) @scale: scaling flag. (1) => max/min are valid.
   * (251) @disp_max: maximum value used for display generation.
   * (252) @disp_min: minimum value used for display generation.
   */
  float lb_f1;
  float pad17;
  float ph0_f1;
  float ph1_f1;
  float max;
  float min;
  float orig_f1;
  float scale;
  float disp_max;
  float disp_min;

  /* (253..255) @pad18 */
  float pad18[3];

  /* (256) @phase2d: 2d plane type code.
   * (257) @x1_f2: second-dimension extraction region origin.
   * (258) @xn_f2: second-dimension extraction region endpoint.
   * (259) @x1_f1: first-dimension extraction region origin.
   * (260) @xn_f1: first-dimension extraction region endpoint.
   * (261) @x1_f3: third-dimension extraction region origin.
   * (262) @xn_f3: third-dimension extraction region endpoint.
   * (263) @x1_f4: fourth-dimension extraction region origin.
   * (264) @xn_f4: fourth-dimension extraction region endpoint.
   */
  float phase2d;
  float x1_f2;
  float xn_f2;
  float x1_f1;
  float xn_f1;
  float x1_f3;
  float xn_f3;
  float x1_f4;
  float xn_f4;

  /* (265..282) @pad19 */
  float pad19[18];

  /* (283) @t_hour: hour of data conversion/collection.
   * (284) @t_min: minute of data conversion/collection.
   * (285) @t_sec: second of data conversion/collection.
   * (286..289) @src: name of source file.
   */
  float t_hour;
  float t_min;
  float t_sec;
  char src[PIPE_HDRSTR_SZ_SRC];

  /* (290..293) @username: user who converted/collected data. */
  char username[PIPE_HDRSTR_SZ_USER];

  /* (294) @d_month: month of data conversion/collection.
   * (295) @d_day: day of data conversion/collection.
   * (296) @d_year: year of data conversion/collection.
   * (297..311) @title: dataset title string.
   */
  float d_month;
  float d_day;
  float d_year;
  char title[PIPE_HDRSTR_SZ_TITLE];

  /* (312..351) @comment: dataset comment string. */
  char comment[PIPE_HDRSTR_SZ_COMMENT];

  /* (352..358) @pad20 */
  float pad20[7];

  /* unused fields appended to spectral data:
   * (359) @block_last
   * (360) @block_cont
   * (361) @block_base
   * (362) @block_peak
   * (363) @block_bmap
   * (364) @block_hist
   * (365) @block_1d
   */
  float block_last;
  float block_cont;
  float block_base;
  float block_peak;
  float block_bmap;
  float block_hist;
  float block_1d;

  /* (366..385) @pad21 */
  float pad21[20];

  /* (386) @tdsz_f2: second-dimension time-domain point count.
   * (387) @tdsz_f1: first-dimension time-domain point count.
   * (388) @tdsz_f3: third-dimension time-domain point count.
   * (389) @tdsz_f4: fourth-dimension time-domain point count.
   */
  float tdsz_f2;
  float tdsz_f1;
  float tdsz_f3;
  float tdsz_f4;

  /* (390..398) @pad22 */
  float pad22[9];

  /* (399) @virgin2d: (0) => data never accessed, header never adjusted.
   * (400) @apodcode_f3: third-dimension window function type.
   * (401..403) @apodq_f3: third-dimension window parameters.
   * (404) @c1_f3: third-dimension first-point scaling modifier.
   * (405) @apodcode_f4: fourth-dimension window function type.
   * (406..408) @apodq_f4: fourth-dimension window parameters.
   * (409) @c1_f4: fourth-dimension first-point scaling modifier.
   */
  float virgin2d;
  float apodcode_f3;
  float apodq_f3[3];
  float c1_f3;
  float apodcode_f4;
  float apodq_f4[3];
  float c1_f4;

  /* (410..412) @pad23 */
  float pad23[3];

  /* (413) @apodcode_f2: second-dimension window function type.
   * (414) @apodcode_f1: first-dimension window function type.
   * (415..417) @apodq_f2: second-dimension window parameters.
   * (418) @c1_f2: second-dimension first-point scaling modifier.
   * (419) @pad24
   * (420..422) @apodq_f1: first-dimension window parameters.
   * (423) @c1_f1: first-dimension first-point scaling modifier.
   */
  float apodcode_f2;
  float apodcode_f1;
  float apodq_f2[3];
  float c1_f2;
  float pad24;
  float apodq_f1[3];
  float c1_f1;

  /* (424..427) @pad25 */
  float pad25[4];

  /* (428) @apod_f1: first-dimension apodization point count. */
  float apod_f1;

  /* (429..436) @pad26 */
  float pad26[8];

  /* (437) @zf_f1: first-dimension zero-fill value.
   * (438) @zf_f3: third-dimension zero-fill value.
   * (439) @zf_f4: fourth-dimension zero-fill value.
   */
  float zf_f1;
  float zf_f3;
  float zf_f4;

  /* (440..441) @pad27 */
  float pad27[2];

  /* (442) @file_count: number of files in the complete datum.
   * (443) @slice_count: number of 1D slices in the stream.
   */
  float file_count;
  float slice_count;

  /* (444..463) @pad28 */
  float pad28[20];

  /* (464..471) @oper: operator name string. */
  char oper[PIPE_HDRSTR_SZ_OPER];

  /* (472..474) @pad29 */
  float pad29[3];

  /* (475) @aqsgn_f1: first-dimension fft sign adjustment.
   * (476) @aqsgn_f3: third-dimension fft sign adjustment.
   * (477) @aqsgn_f4: fourth-dimension fft sign adjustment.
   */
  float aqsgn_f1;
  float aqsgn_f3;
  float aqsgn_f4;

  /* (478..479) @pad30 */
  float pad30[2];

  /* (480) @offppm_f2: second-dimension additional alignment offset.
   * (481) @offppm_f1: first-dimension additional alignment offset.
   * (482) @offppm_f3: third-dimension additional alignment offset.
   * (483) @offppm_f4: fourth-dimension additional alignment offset.
   */
  float offppm_f2;
  float offppm_f1;
  float offppm_f3;
  float offppm_f4;

  /* (484..511) @pad_end */
  float pad_end[28];
};

/* pipe_read_header(): read the contents of a header from a pipe-format file.
 * @fname: the input filename.
 * @endianness: byte order result pointer.
 * @hdr: header result pointer.
 */
int pipe_read_header (const char *fname,
                      enum byteorder *endianness,
                      struct pipe_header *hdr) {
  /* declare a few required variables:
   */
  unsigned int n_bytes, n_words;
  uint8_t *bytes;

  /* read in the file header bytes. */
  n_bytes = sizeof(struct pipe_header);
  n_words = n_bytes / sizeof(float);
  bytes = bytes_read_block(fname, 0, n_bytes);

  /* check that the bytes were read successfully. */
  if (!bytes)
    throw("failed to read header from '%s'", fname);

  /* copy the header bytes onto the header structure. */
  memcpy(hdr, bytes, n_bytes);

  /* check if the floating point order value is correct. if not, the bytes
   * likely need to be swapped.
   */
  if (hdr->order != 0.0 && hdr->order != (float) PIPE_MAGIC) {
    /* swap the bytes of each word. */
    bytes_swap(bytes, n_words, sizeof(float));

    /* re-copy the bytes onto the structure. */
    memcpy(hdr, bytes, n_bytes);

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

/* pipe_interlace(): pipe-format files do not store their individual traces
 * as interlaced R,I,R,I..., but instead interlaces entire real traces and
 * imaginary traces. this fixes pipe's mistake so hx_array_complexify()
 * can operate on pipe data just like any other data type.
 * @x: pointer to the array to interlace.
 * @n: number of *complex* points per trace.
 */
int pipe_interlace (hx_array *x, unsigned int n) {
  /* declare a few required variables:
   * @i: trace loop counter.
   * @j: point loop counter.
   * @ntraces: number of complex traces.
   * @kre: real trace start index.
   * @kim: imag trace start index.
   * @xtmp: temporary duplicate coefficient array.
   */
  unsigned int i, j, ntraces, isrc_re, isrc_im, idest_re, idest_im;
  real *xtmp;

  /* check that the trace size evenly divides the array elements. */
  if (x->sz[0] % (2 * n))
    throw("trace size %u does not evenly divide array (%d)", n, x->sz[0]);

  /* compute the number of traces in the array. */
  ntraces = x->sz[0] / (2 * n);

  /* allocate a temporary duplicate array of the array elements. */
  xtmp = (real*) malloc(x->len * sizeof(real));
  if (!xtmp)
    throw("failed to allocate temporary array of %d reals", x->len);

  /* copy the data over into the duplicate array. */
  memcpy(xtmp, x->x, x->len * sizeof(real));

  /* loop over each trace in the linear real array. */
  for (i = 0; i < ntraces; i++) {
    /* compute the real and imaginary trace starting indices. */
    isrc_re = (2 * i) * n;
    isrc_im = (2 * i + 1) * n;

    /* loop over the points of the trace. */
    for (j = 0; j < n; j++) {
      /* compute the output indices. */
      idest_re = i * (2 * n) + 2 * j;
      idest_im = i * (2 * n) + 2 * j + 1;

      /* shuffle the points around. */
      x->x[idest_re] = xtmp[isrc_re + j];
      x->x[idest_im] = xtmp[isrc_im + j];
    }
  }

  /* free the duplicate array. */
  free(xtmp);

  /* return success. */
  return 1;
}

/* pipe_guess(): check whether a file contains pipe-format data.
 * @fname: the input filename.
 */
int pipe_guess (const char *fname) {
  /* declare a few required variables:
   * @fh: input file handle.
   * @wd: first words of input file.
   */
  float wd[3];
  FILE *fh;

  /* open the input file. */
  fh = fopen(fname, "rb");

  /* check that the file was opened. */
  if (!fh)
    throw("failed to open '%s'", fname);

  /* read the first three words. */
  if (!fread(&wd, sizeof(float), 3, fh))
    throw("failed to read magic numbers");

  /* close the input file. */
  fclose(fh);

  /* check the magic word, without swapping. */
  if (wd[2] == (float) PIPE_MAGIC)
    return 1;

  /* swap the word and check again. */
  bytes_swap_u32((uint32_t*) &wd[2]);
  if (wd[2] == (float) PIPE_MAGIC)
    return 1;

  /* no match. */
  return 0;
}

/* pipe_decode(): read pipe parameters into a datum structure.
 * @D: pointer to the destination datum structure.
 * @fname: the input filename.
 */
int pipe_decode (datum *D, const char *fname) {
  /* declare variables required to determine byte ordering:
   * @endian: the byte ordering of the data file.
   * @hdr: the pipe file header structure.
   */
  enum byteorder endian = BYTES_ENDIAN_AUTO;
  struct pipe_header hdr;

  /* declare variables required to traverse dimensions:
   * @d: dimension loop counter.
   * @ord: dimension index array.
   */
  unsigned int d;
  int ord[PIPE_MAXDIM];

  /* declare variables required to parse date information:
   * @ts: calendar time structure to hold header values.
   */
  struct tm ts;

  /* read the header information from the data file. */
  if (!pipe_read_header(fname, &endian, &hdr))
    throw("failed to read header of '%s'", fname);

  /* check if time/date information exists in the header fields. */
  if (hdr.d_year && hdr.d_month) {
    /* parse the time header fields. */
    ts.tm_hour = (int) hdr.t_hour;
    ts.tm_min = (int) hdr.t_min;
    ts.tm_sec = (int) hdr.t_sec;

    /* parse the date header fields. */
    ts.tm_mon = (int) hdr.d_month - 1;
    ts.tm_mday = (int) hdr.d_day;
    ts.tm_year = (int) hdr.d_year - 1900;

    /* convert the date and time into an epoch offset. */
    D->epoch = mktime(&ts);
  }

  /* check the dimensionality. */
  if ((int) hdr.ndims < 1 ||
      (int) hdr.ndims > PIPE_MAXDIM)
    throw("invalid dimensionality %.0f", hdr.ndims);

  /* initially set the number of dimensions to the maximum allowed, because
   * pipe arranges its dimension information in a really screwy way.
   */
  if (!datum_dims_realloc(D, PIPE_MAXDIM))
    throw("failed to allocate dimension array");

  /* store the dimension ordering array values. */
  for (d = 0; d < PIPE_MAXDIM; d++)
    ord[d] = hdr.dimorder[d] - 1;

  /* store the complex states from the quadrature flags. */
  D->dims[ord[0]].cx = ((int) hdr.quad_f1 != PIPE_QUAD_REAL);
  D->dims[ord[1]].cx = ((int) hdr.quad_f2 != PIPE_QUAD_REAL);
  D->dims[ord[2]].cx = ((int) hdr.quad_f3 != PIPE_QUAD_REAL);
  D->dims[ord[3]].cx = ((int) hdr.quad_f4 != PIPE_QUAD_REAL);

  /* store the alternation states from the sign flags. */
  D->dims[ord[0]].alt = ((int) hdr.aqsgn_f1 & PIPE_AQSGN_ALT ? 1 : 0);
  D->dims[ord[1]].alt = ((int) hdr.aqsgn_f2 & PIPE_AQSGN_ALT ? 1 : 0);
  D->dims[ord[2]].alt = ((int) hdr.aqsgn_f3 & PIPE_AQSGN_ALT ? 1 : 0);
  D->dims[ord[3]].alt = ((int) hdr.aqsgn_f4 & PIPE_AQSGN_ALT ? 1 : 0);

  /* store the negation states from the sign flags. */
  D->dims[ord[0]].neg = ((int) hdr.aqsgn_f1 & PIPE_AQSGN_NEG ? 1 : 0);
  D->dims[ord[1]].neg = ((int) hdr.aqsgn_f2 & PIPE_AQSGN_NEG ? 1 : 0);
  D->dims[ord[2]].neg = ((int) hdr.aqsgn_f3 & PIPE_AQSGN_NEG ? 1 : 0);
  D->dims[ord[3]].neg = ((int) hdr.aqsgn_f4 & PIPE_AQSGN_NEG ? 1 : 0);

  /* store the gradient-enhanced states from the quadrature flags. */
  D->dims[ord[0]].genh = ((int) hdr.quad_f1 == PIPE_QUAD_GE);
  D->dims[ord[1]].genh = ((int) hdr.quad_f2 == PIPE_QUAD_GE);
  D->dims[ord[2]].genh = ((int) hdr.quad_f3 == PIPE_QUAD_GE);
  D->dims[ord[3]].genh = ((int) hdr.quad_f4 == PIPE_QUAD_GE);

  /* store the fourier-transform flags. */
  D->dims[ord[0]].ft = (unsigned int) hdr.ftflag_f1;
  D->dims[ord[1]].ft = (unsigned int) hdr.ftflag_f2;
  D->dims[ord[2]].ft = (unsigned int) hdr.ftflag_f3;
  D->dims[ord[3]].ft = (unsigned int) hdr.ftflag_f4;

  /* store the nucleus strings. */
  strncpy(D->dims[ord[0]].nuc, hdr.label_f1, PIPE_HDRSTR_SZ_LABEL);
  strncpy(D->dims[ord[1]].nuc, hdr.label_f2, PIPE_HDRSTR_SZ_LABEL);
  strncpy(D->dims[ord[2]].nuc, hdr.label_f3, PIPE_HDRSTR_SZ_LABEL);
  strncpy(D->dims[ord[3]].nuc, hdr.label_f4, PIPE_HDRSTR_SZ_LABEL);

  /* null-terminate the nucleus string. */
  for (d = 0; d < PIPE_MAXDIM; d++)
    D->dims[d].nuc[7] = '\0';

  /* store the size parameters. */
  D->dims[ord[0]].td = D->dims[ord[0]].tdunif = hdr.tdsz_f1;
  D->dims[ord[1]].td = D->dims[ord[1]].tdunif = hdr.tdsz_f2;
  D->dims[ord[2]].td = D->dims[ord[2]].tdunif = hdr.tdsz_f3;
  D->dims[ord[3]].td = D->dims[ord[3]].tdunif = hdr.tdsz_f4;

  /* store the f1 current size parameter. */
  D->dims[ord[0]].sz = hdr.x1_f1 && hdr.xn_f1 ? hdr.xn_f1 - hdr.x1_f1 + 1 :
                       D->dims[ord[0]].ft ? hdr.ftsz_f1 : hdr.apod_f1;

  /* store the f2 current size parameter. */
  D->dims[ord[1]].sz = hdr.x1_f2 && hdr.xn_f2 ? hdr.xn_f2 - hdr.x1_f2 + 1 :
                       D->dims[ord[1]].ft ? hdr.ftsz_f2 : hdr.apod_f2;

  /* store the f3 current size parameter. */
  D->dims[ord[2]].sz = hdr.x1_f3 && hdr.xn_f3 ? hdr.xn_f3 - hdr.x1_f3 + 1 :
                       D->dims[ord[2]].ft ? hdr.ftsz_f3 : hdr.apod_f3;

  /* store the f4 current size parameter. */
  D->dims[ord[3]].sz = hdr.x1_f4 && hdr.xn_f4 ? hdr.xn_f4 - hdr.x1_f4 + 1 :
                       D->dims[ord[3]].ft ? hdr.ftsz_f4 : hdr.apod_f4;

  /* store the spectral width parameters. */
  D->dims[ord[0]].width = hdr.sw_f1;
  D->dims[ord[1]].width = hdr.sw_f2;
  D->dims[ord[2]].width = hdr.sw_f3;
  D->dims[ord[3]].width = hdr.sw_f4;

  /* store the carrier frequency parameters. */
  D->dims[ord[0]].carrier = hdr.obs_f1;
  D->dims[ord[1]].carrier = hdr.obs_f2;
  D->dims[ord[2]].carrier = hdr.obs_f3;
  D->dims[ord[3]].carrier = hdr.obs_f4;

  /* store the spectral offset parameters. */
  D->dims[ord[0]].offset = hdr.car_f1 * hdr.obs_f1;
  D->dims[ord[1]].offset = hdr.car_f2 * hdr.obs_f2;
  D->dims[ord[2]].offset = hdr.car_f3 * hdr.obs_f3;
  D->dims[ord[3]].offset = hdr.car_f4 * hdr.obs_f4;

  /* set the true dimension count and reallocate the dimension array. */
  if (!datum_dims_realloc(D, hdr.ndims))
    throw("failed to reallocate dimension array");

  /* store the filename string. */
  D->fname = (char*) malloc((strlen(fname) + 1) * sizeof(char));
  if (D->fname)
    strcpy(D->fname, fname);

  /* store the datum type. */
  D->type = DATUM_TYPE_PIPE;
  D->endian = endian;

  /* return success. */
  return 1;
}

/* pipe_fwrite_dim(): write the contents of a datum core array into an
 * output file stream of pipe-format coefficients. recursively drop from
 * cubes, to planes, to traces, and finally to points.
 * @D: pointer to the source datum structure.
 * @dim: current datum dimension in the recursion.
 * @n0: current coefficient offset in the recursion.
 * @arr: index array for looping through the datum core array.
 * @fh: the output file handle.
 */
int pipe_fwrite_dim (datum *D, unsigned int dim, int n0, int *arr, FILE *fh) {
  /* declare a few required variables:
   * @d: current array algebraic dimension index.
   * @k: current array topological dimension index.
   * @n: current array coefficient imaginary offset.
   * @num: size along current array dimension.
   * @idx: array linear scalar index.
   * @f: float output value.
   */
  int d, k, n, num, idx;
  float f;

  /* determine the array dimension indices. */
  d = D->dims[dim].d;
  k = D->dims[dim].k;

  /* determine the coefficient index. */
  n = (d == DATUM_DIM_INVALID ? 0 : 1 << d);

  /* compute the number of points in the dimension. */
  num = D->array.sz[k];

  /* check if we've reached the lowest dimension. */
  if (k == 0) {
    /* store the real points of the current trace. */
    for (arr[k] = 0; arr[k] < num; arr[k]++) {
      /* pack the linear index and get the coefficient. */
      hx_array_index_pack(D->array.k, D->array.sz, arr, &idx);
      f = (float) D->array.x[D->array.n * idx + n0];

      /* write the coefficient. */
      if (fwrite(&f, sizeof(float), 1, fh) != 1)
        return 0;
    }

    /* check if the trace is complex. */
    if (D->dims[dim].cx) {
      /* store the imaginary points of the current trace. */
      for (arr[k] = 0; arr[k] < num; arr[k]++) {
        /* pack the linear index and get the coefficient. */
        hx_array_index_pack(D->array.k, D->array.sz, arr, &idx);
        f = (float) D->array.x[D->array.n * idx + n0 + n];

        /* write the coefficient. */
        if (fwrite(&f, sizeof(float), 1, fh) != 1)
          return 0;
      }
    }
  }
  else {
    /* recurse into the lower dimensions. */
    for (arr[k] = 0; arr[k] < num; arr[k]++) {
      /* recurse the real component. */
      if (!pipe_fwrite_dim(D, dim - 1, n0, arr, fh))
        return 0;

      /* recurse the imaginary component. */
      if (D->dims[dim].cx && !pipe_fwrite_dim(D, dim - 1, n0 + n, arr, fh))
        return 0;
    }
  }

  /* return success. */
  return 1;
}

/* pipe_encode(): write a datum structure in pipe-format to a file.
 * @D: pointer to the source structure.
 * @fname: the output filename.
 */
int pipe_encode (datum *D, const char *fname) {
  /* declare variables required to output pipe-format files:
   * @ord: dimension ordering array.
   * @arr: array index for iteration during output.
   * @hdr: the pipe file header structure.
   * @nhdr: number of floats in @fhdr.
   * @ts: calendar time structure.
   * @fh: the output file handle.
   */
  int ord[PIPE_MAXDIM], arr[PIPE_MAXDIM];
  struct pipe_header hdr;
  int i, n, nhdr;
  struct tm *ts;
  FILE *fh;

  /* check that the datum will fit in a pipe-format file. */
  if (D->nd > PIPE_MAXDIM)
    throw("datum contains too many dimensions for pipe format");

  /* initialize the values in the header. */
  memset(&hdr, 0, sizeof(struct pipe_header));

  /* set the magic numbers. */
  hdr.magic = 0.0;
  hdr.format = (float) 0xeeeeeeee;
  hdr.order = (float) PIPE_MAGIC;

  /* compute the calendar date structure fields. */
  ts = gmtime(&D->epoch);

  /* store the header date fields. */
  hdr.d_year = (float) ts->tm_year + 1900;
  hdr.d_month = (float) ts->tm_mon + 1.0;
  hdr.d_day = (float) ts->tm_mday;

  /* store the header time fields. */
  hdr.t_hour = (float) ts->tm_hour;
  hdr.t_min = (float) ts->tm_min;
  hdr.t_sec = (float) ts->tm_sec;

  /* set the number of dimensions and the plane mode. */
  hdr.ndims = (float) D->nd;
  hdr.phase2d = PIPE_PHASE2D_STATES;

  /* set the master quadrature flag. */
  hdr.quad = (float) PIPE_QUAD_REAL;
  for (i = 0; i < D->nd; i++) {
    /* having even one complex dimension sets the complex flag. */
    if (D->dims[i].cx || (int) hdr.quad == PIPE_QUAD_COMPLEX)
      hdr.quad = PIPE_QUAD_COMPLEX;
  }

  /* set the dimension order (YXZA). */
  for (i = 0; i < PIPE_MAXDIM; i++)
    hdr.dimorder[i] = (float) (ord[i] = i + 1);

  /* if necessary, write the first-dimension header information. */
  if (D->nd >= 1) {
    /* set quadrature information. */
    hdr.quad_f1 = (D->dims[0].cx ? PIPE_QUAD_COMPLEX : PIPE_QUAD_REAL);
    hdr.aqsgn_f1 = PIPE_AQSGN_NONE;
    hdr.ftflag_f1 = (float) D->dims[0].ft;

    /* set nucleus information. */
    strncpy(hdr.label_f1, D->dims[0].nuc, PIPE_HDRSTR_SZ_LABEL);

    /* set size information. */
    hdr.tdsz_f1 = (float) D->dims[0].td;
    hdr.ftsz_f1 = (float) D->dims[0].sz;
    hdr.apod_f1 = (float) D->dims[0].sz;
    hdr.sz = (float) D->dims[0].sz;

    /* set spectral width, offset and carrier. */
    hdr.sw_f1 = (float) D->dims[0].width;
    hdr.obs_f1 = (float) D->dims[0].carrier;
    hdr.orig_f1 = (float) D->dims[0].offset;
  }

  /* if necessary, write the second-dimension header information. */
  if (D->nd >= 2) {
    /* set quadrature information. */
    hdr.quad_f2 = (D->dims[1].cx ? PIPE_QUAD_COMPLEX : PIPE_QUAD_REAL);
    hdr.aqsgn_f2 = PIPE_AQSGN_NONE;
    hdr.ftflag_f2 = (float) D->dims[1].ft;

    /* set nucleus information. */
    strncpy(hdr.label_f2, D->dims[1].nuc, PIPE_HDRSTR_SZ_LABEL);

    /* set size information. */
    hdr.tdsz_f2 = (float) D->dims[1].td;
    hdr.ftsz_f2 = (float) D->dims[1].sz;
    hdr.apod_f2 = (float) D->dims[1].sz;
    hdr.specnum = (float) D->dims[1].sz;
    if (D->dims[1].cx) hdr.specnum *= 2.0;

    /* set spectral width, offset and carrier. */
    hdr.sw_f2 = (float) D->dims[1].width;
    hdr.obs_f2 = (float) D->dims[1].carrier;
    hdr.orig_f2 = (float) D->dims[1].offset -  hdr.sw_f2 / 2.0;
  }

  /* if necessary, write the third-dimension header information. */
  if (D->nd >= 3) {
    /* set the pipe flag. */
    hdr.pipe = 1;

    /* set quadrature information. */
    hdr.quad_f3 = (D->dims[2].cx ? PIPE_QUAD_COMPLEX : PIPE_QUAD_REAL);
    hdr.aqsgn_f3 = PIPE_AQSGN_NONE;
    hdr.ftflag_f3 = (float) D->dims[2].ft;

    /* set nucleus information. */
    strncpy(hdr.label_f3, D->dims[2].nuc, PIPE_HDRSTR_SZ_LABEL);

    /* set size information. */
    hdr.tdsz_f3 = (float) D->dims[2].td;
    hdr.ftsz_f3 = (float) D->dims[2].sz;
    hdr.apod_f3 = (float) D->dims[2].sz;
    hdr.size_f3 = (float) D->dims[2].sz;
    if (D->dims[2].cx) hdr.size_f3 *= 2.0;

    /* set spectral width, offset and carrier. */
    hdr.sw_f3 = (float) D->dims[2].width;
    hdr.obs_f3 = (float) D->dims[2].carrier;
    hdr.orig_f3 = (float) D->dims[2].offset -  hdr.sw_f3 / 2.0;
  }

  /* if necessary, write the fourth-dimension header information. */
  if (D->nd >= 4) {
    /* set quadrature information. */
    hdr.quad_f4 = (D->dims[3].cx ? PIPE_QUAD_COMPLEX : PIPE_QUAD_REAL);
    hdr.aqsgn_f4 = PIPE_AQSGN_NONE;
    hdr.ftflag_f4 = (float) D->dims[3].ft;

    /* set nucleus information. */
    strncpy(hdr.label_f4, D->dims[3].nuc, PIPE_HDRSTR_SZ_LABEL);

    /* set size information. */
    hdr.tdsz_f4 = (float) D->dims[3].td;
    hdr.ftsz_f4 = (float) D->dims[3].sz;
    hdr.apod_f4 = (float) D->dims[3].sz;
    hdr.size_f4 = (float) D->dims[3].sz;
    if (D->dims[3].cx) hdr.size_f4 *= 2.0;

    /* set spectral width, offset and carrier. */
    hdr.sw_f4 = (float) D->dims[3].width;
    hdr.obs_f4 = (float) D->dims[3].carrier;
    hdr.orig_f4 = (float) D->dims[3].offset -  hdr.sw_f4 / 2.0;
  }

  /* compute the size of the header buffer. */
  nhdr = sizeof(struct pipe_header) / sizeof(float);

  /* initialize the data buffer index. */
  arr[0] = arr[1] = arr[2] = arr[3] = 0;
  n = 0;

  /* open the output file. */
  if (fname)
    fh = fopen(fname, "wb");
  else
    fh = stdout;

  /* check that the output file was opened. */
  if (!fh)
    throw("failed to open '%s'", fname);

  /* write the header to the output stream. */
  if (fwrite(&hdr, sizeof(float), nhdr, fh) != nhdr)
    throw("failed to write %d header values", nhdr);

  /* write coefficients out from the datum core array. */
  if (!pipe_fwrite_dim(D, D->nd - 1, n, arr, fh))
    throw("failed to convert core array to pipe format");

  /* close the output file, if necessary. */
  if (fname)
    fclose(fh);

  /* return success. */
  return 1;
}

/* pipe_array(): read a pipe-format data file into a datum array.
 * @D: pointer to the destination datum structure.
 */
int pipe_array (datum *D) {
  /* declare a few required variables:
   * @endian: byte ordering of the input file.
   * @i: dimension loop counter.
   * @n_words: number of data words in the file.
   * @n_expected: expected file size.
   * @n_actual: actual file size.
   * @offset: header byte offset.
   * @hdr: pipe header structure.
   * @fh: input file handle.
   */
  enum byteorder endian = BYTES_ENDIAN_AUTO;
  unsigned int d, n_words, n_expected, n_actual, offset;
  struct pipe_header hdr;
  FILE *fh;

  /* check that the input filename is valid. */
  if (D->fname == NULL)
    throw("invalid input filename");

  /* read the header information from the data file. */
  if (!pipe_read_header(D->fname, &endian, &hdr))
    throw("failed to read header of '%s'", D->fname);

  /* compute the byte offset from which to begin reading point data. */
  offset = sizeof(struct pipe_header);

  /* compute the actual number of bytes in the file. */
  n_actual = bytes_size(D->fname) - offset;
  n_words = n_actual / sizeof(float);

  /* calculate the expected number of bytes in the file. */
  for (d = 0, n_expected = sizeof(float); d < D->nd; d++)
    n_expected *= (D->dims[d].cx ? 2 : 1) * D->dims[d].sz;

  /* check that the byte counts match. */
  if (n_expected != n_actual)
    throw("data size mismatch (expected %u, read %u)", n_expected, n_actual);

  /* open the input file for reading. */
  fh = fopen(D->fname, "rb");

  /* check that the file was opened. */
  if (!fh)
    throw("failed to open '%s'", D->fname);

  /* read data from the file into the output array. */
  if (!hx_array_fread_raw(fh, &D->array, endian, sizeof(float), 1,
                          offset, 0, 1, n_words, 0))
    throw("failed to read raw data from '%s'", D->fname);

  /* close the input file. */
  fclose(fh);

  /* interlace the complex direct-dimension traces, if necessary. */
  if (D->dims[0].cx && !pipe_interlace(&D->array, D->dims[0].sz))
    throw("failed to interlace complex traces");

  /* return success. */
  return 1;
}

