
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
#ifndef __HXND_NMR_PIPE_H__
#define __HXND_NMR_PIPE_H__

/* include the n-dimensional math header. */
#include <hxnd/hx.h>

/* include the byte-level data and nmr datum headers. */
#include <hxnd/nmr-datum.h>
#include <hxnd/bytes.h>

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
#define PIPE_2D_STATES  2

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

/* function declarations: */

int pipe_check_magic (const char *fname);

int pipe_read_header (const char *fname, enum byteorder *endianness,
                      struct pipe_header *hdr);

int pipe_read (const char *fname, unsigned int n, hx_array *x);

int pipe_interlace (hx_array *x, unsigned int n);

int pipe_fill_datum (const char *fname, datum *D);

int pipe_fwrite_datum (datum *D, FILE *fh);

#endif /* __HXND_NMR_PIPE_H__ */

