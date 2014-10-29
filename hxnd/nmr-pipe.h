
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

/* pipe_header: properly padded header of data contained in a pipe-format file.
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
   */
  float ndims;

  /* (10) @obs_f3: third-dimension observe frequency.
   * (11) @sw_f3: third-dimension spectral witdh.
   * (12) @orig_f3: third-dimension spectral origin.
   * (13) @ftflag_f3: third-dimension fft status.
   */
  float obs_f3;
  float sw_f3;
  float orig_f3;
  float ftflag_f3;

  /* (14) @plane_loc: FIXME
   */
  float plane_loc;

  /* (15) @size_f3: third-dimension size.
   */
  float size_f3;

  /* (16..17) @label_f2: second-dimension label.
   * (18..19) @label_f1: first-dimension label.
   * (20..21) @label_f3: third-dimension label.
   * (22..23) @label_f4: fourth-dimension label.
   */
  float label_f2[2];
  float label_f1[2];
  float label_f3[2];
  float label_f4[2];

  /* (24) @dimorder_f1:
   * (25) @dimorder_f2:
   * (26) @dimorder_f3:
   * (27) @dimorder_f4:
   */
  float dimorder_f1;
  float dimorder_f2;
  float dimorder_f3;
  float dimorder_f4;

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

  /* (50) @apod_f3: third-dimension apodization flag.
   * (51) @quad_f3: third-dimension quadrature flag.
   */
  float apod_f3;
  float quad_f3;

  /* (52) @pad2 */
  float pad2;

  /* (53) @apod_f4: fourth-dimension apodization flag.
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
   * (64) @aqsgn_f2: second-dimension acquisition sign.
   */
  float pipe;
  float units_f3;
  float units_f4;
  float ph0_f3;
  float ph1_f3;
  float ph0_f4;
  float ph1_f4;
  float aqsgn_f2;

  /* (65) @partition: FIXME
   */
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
   * (75) @pipecount: FIXME
   * (76) @pad3
   */
  float user[5];
  float pipecount;
  float pad3;

  /* (77) @first_plane: FIXME
   * (78) @last_plane: FIXME
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

  /* (95) @apod_f2: second-dimension apodization flag.
   * (96) @ftsz_f2: second-dimension fft points count.
   * (97) @realsz: FIXME
   * (98) @ftsz_f1: first-dimension fft points count.
   * (99) @sz: FIXME
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

  /* (119) @obs_f2: second-dimension observe frequency.
   */
  float obs_f2;

  /* (120..134) @pad8 */
  float pad8[15];

  /* (135) @mcflag: FIXME
   */
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

  /* (157) @temperature: sample temperature. */
  float temperature;

  /* (158..179) @pad11 */
  float pad11[22];

  /* (180) @rank: FIXME
   */
  float rank;

  /* (181..198) @pad12 */
  float pad12[18];

  /* (199) @tau: FIXME
   * (200) @ftsz_f3: third-dimension fft points count.
   * (201) @ftsz_f4: fourth-dimension fft points count.
   */
  float tau;
  float ftsz_f3;
  float ftsz_f4;

  /* (202..217) @pad13 */
  float pad13[16];

  /* (218) @obs_f1: first-dimension observe frequency.
   * (219) @specnum: FIXME
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
   * (247) @max: FIXME
   * (248) @min: FIXME
   * (249) @orig_f1: first-dimension spectral origin.
   * (250) @scale: scaling flag.
   * (251) @disp_max: FIXME
   * (252) @disp_min: FIXME
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

  /* (256) @phase2d: FIXME
   * (257) @x1_f2: FIXME
   * (258) @xn_f2: FIXME
   * (259) @x1_f1: FIXME
   * (260) @xn_f1: FIXME
   * (261) @x1_f3: FIXME
   * (262) @xn_f3: FIXME
   * (263) @x1_f4: FIXME
   * (264) @xn_f4: FIXME
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

  /* (283) @t_hour: FIXME
   * (284) @t_min: FIXME
   * (285) @t_sec: FIXME
   * (286) @src: FIXME
   */
  float t_hour;
  float t_min;
  float t_sec;
  float src;

  /* (287..289) @pad20 */
  float pad20[3];

  /* (290) @username: FIXME */
  float username;

  /* (291..293) @pad21 */
  float pad21[3];

  /* (294) @d_month: FIXME
   * (295) @d_day: FIXME
   * (296) @d_year: FIXME
   * (297) @title: FIXME
   */
  float d_month;
  float d_day;
  float d_year;
  float title;

  /* (298..311) @pad22 */
  float pad22[14];

  /* (312) @comment: user comment. */
  float comment;

  /* (313..358) @pad23 */
  float pad23[46];

  /* (359) @block_last: FIXME
   * (360) @block_cont: FIXME
   * (361) @block_base: FIXME
   * (362) @block_peak: FIXME
   * (363) @block_bmap: FIXME
   * (364) @block_hist: FIXME
   * (365) @block_1d: FIXME
   */
  float block_last;
  float block_cont;
  float block_base;
  float block_peak;
  float block_bmap;
  float block_hist;
  float block_1d;

  /* (366..385) @pad24 */
  float pad24[20];

  /* (386) @tdsz_f2: second-dimension time-domain point count.
   * (387) @tdsz_f1: first-dimension time-domain point count.
   * (388) @tdsz_f3: third-dimension time-domain point count.
   * (389) @tdsz_f4: fourth-dimension time-domain point count.
   */
  float tdsz_f2;
  float tdsz_f1;
  float tdsz_f3;
  float tdsz_f4;

  /* (390..398) @pad25 */
  float pad25[9];

  /* (399) @virgin2d: FIXME
   * (400) @apodcode_f3: FIXME
   * (401) @apodq1_f3: FIXME
   * (402) @apodq2_f3: FIXME
   * (403) @apodq3_f3: FIXME
   * (404) @c1_f3: FIXME
   * (405) @apodcode_f4: FIXME
   * (406) @apodq1_f4: FIXME
   * (407) @apodq2_f4: FIXME
   * (408) @apodq3_f4: FIXME
   * (409) @c1_f4: FIXME
   */
  float virgin2d;
  float apodcode_f3;
  float apodq1_f3;
  float apodq2_f3;
  float apodq3_f3;
  float c1_f3;
  float apodcode_f4;
  float apodq1_f4;
  float apodq2_f4;
  float apodq3_f4;
  float c1_f4;

  /* (410..412) @pad26 */
  float pad26[3];

  /* (413) @apodcode_f2: FIXME
   * (414) @apodcode_f1: FIXME
   * (415) @apodq1_f2: FIXME
   * (416) @apodq2_f2: FIXME
   * (417) @apodq3_f2: FIXME
   * (418) @c1_f2: FIXME
   * (419) @pad27
   * (420) @apodq1_f1: FIXME
   * (421) @apodq2_f1: FIXME
   * (422) @apodq3_f1: FIXME
   * (423) @c1_f1: FIXME
   */
  float apodcode_f2;
  float apodcode_f1;
  float apodq1_f2;
  float apodq2_f2;
  float apodq3_f2;
  float c1_f2;
  float pad27;
  float apodq1_f1;
  float apodq2_f1;
  float apodq3_f1;
  float c1_f1;

  /* (424..427) @pad28 */
  float pad28[4];

  /* (428) @apod_f1: first-dimension apodization flag. */
  float apod_f1;

  /* (429..436) @pad29 */
  float pad29[8];

  /* (437) @zf_f1: first-dimension zero-fill value.
   * (438) @zf_f3: third-dimension zero-fill value.
   * (439) @zf_f4: fourth-dimension zero-fill value.
   */
  float zf_f1;
  float zf_f3;
  float zf_f4;

  /* (440..441) @pad30 */
  float pad30[2];

  /* FIXME */

  /* (484..511) @pad_end */
  float pad_end[28];
};

/* function declarations: */

int pipe_read_header (const char *fname, enum byteorder *endianness,
                      struct pipe_header *hdr);

int pipe_fill_datum (const char *fname, datum *D);

#endif /* __HXND_NMR_PIPE_H__ */

