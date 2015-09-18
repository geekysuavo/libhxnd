
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
#ifndef __HXND_FN_ARGS_H__
#define __HXND_FN_ARGS_H__

/* include the processing function header. */
#include <hxnd/fn.h>

static fn_arg fn_args_add[] = {
  { "const",    { .f = 0.0  }, 0, FN_VALTYPE_FLOAT },
  { "file",     { .s = NULL }, 0, FN_VALTYPE_STRING },
  { "scale",    { .f = 1.0  }, 0, FN_VALTYPE_FLOAT },
  { "subtract", { .b = 0    }, 0, FN_VALTYPE_BOOL },
  { NULL,       {},            0, FN_VALTYPE_UNKNOWN }
};

static fn_arg fn_args_baseline[] = {
  { "smooth", { .f = 1.0 }, 0, FN_VALTYPE_FLOAT },
  { NULL,     {},           0, FN_VALTYPE_UNKNOWN }
};

static fn_arg fn_args_crop[] = {
  { "from", { .f = 0.0 }, 0, FN_VALTYPE_FLOAT },
  { "to",   { .f = 1.0 }, 0, FN_VALTYPE_FLOAT },
  { "ppm",  { .b = 0   }, 0, FN_VALTYPE_BOOL },
  { "hz",   { .b = 0   }, 0, FN_VALTYPE_BOOL },
  { NULL,   {},           0, FN_VALTYPE_UNKNOWN }
};

static fn_arg fn_args_cut[] = {
  { "trace", { .iv = NULL }, 0, FN_VALTYPE_INTS },
  { "plane", { .iv = NULL }, 0, FN_VALTYPE_INTS },
  { NULL,    {},             0, FN_VALTYPE_UNKNOWN }
};

static fn_arg fn_args_ffm[] = {
  { "func",  { .s = NULL }, 0, FN_VALTYPE_STRING },
  { "iters", { .i = 1000 }, 0, FN_VALTYPE_INT },
  { NULL,    {},            0, FN_VALTYPE_UNKNOWN }
};

static fn_arg fn_args_fft[] = {
  { "alternate", { .b = 0 }, 0, FN_VALTYPE_BOOL },
  { "negate",    { .b = 0 }, 0, FN_VALTYPE_BOOL },
  { "inverse",   { .b = 0 }, 0, FN_VALTYPE_BOOL },
  { NULL,        {},         0, FN_VALTYPE_UNKNOWN }
};

static fn_arg fn_args_filter[] = {
  { "order", { .i = 32       }, 0, FN_VALTYPE_INT },
  { "lo",    { .f = INFINITY }, 0, FN_VALTYPE_FLOAT },
  { "hi",    { .f = INFINITY }, 0, FN_VALTYPE_FLOAT },
  { "ppm",   { .b = 0        }, 0, FN_VALTYPE_BOOL },
  { "hz",    { .b = 0        }, 0, FN_VALTYPE_BOOL },
  { NULL,    {},                0, FN_VALTYPE_UNKNOWN }
};

static fn_arg fn_args_irls[] = {
  { "norm",  { .f = 1.0 }, 0, FN_VALTYPE_FLOAT },
  { "iters", { .i = 10  }, 0, FN_VALTYPE_INT },
  { NULL,    {},           0, FN_VALTYPE_UNKNOWN }
};

static fn_arg fn_args_ist[] = {
  { "thresh", { .f = 0.9 }, 0, FN_VALTYPE_FLOAT },
  { "iters",  { .i = 200 }, 0, FN_VALTYPE_INT },
  { NULL,     {},           0, FN_VALTYPE_UNKNOWN }
};

static fn_arg fn_args_multiply[] = {
  { "first",  { .f = 0.0 }, 0, FN_VALTYPE_FLOAT },
  { "factor", { .f = 1.0 }, 0, FN_VALTYPE_FLOAT },
  { "invert", { .b = 0   }, 0, FN_VALTYPE_BOOL },
  { NULL,     {},           0, FN_VALTYPE_UNKNOWN }
};

static fn_arg fn_args_phase[] = {
  { "ph0",     { .f = 0.0 }, 0, FN_VALTYPE_FLOAT },
  { "ph1",     { .f = 0.0 }, 0, FN_VALTYPE_FLOAT },
  { "pivot",   { .f = 0.5 }, 0, FN_VALTYPE_FLOAT },
  { "ppm",     { .b = 0   }, 0, FN_VALTYPE_BOOL },
  { "hz",      { .b = 0   }, 0, FN_VALTYPE_BOOL },
  { "inverse", { .b = 0   }, 0, FN_VALTYPE_BOOL },
  { NULL,      {},           0, FN_VALTYPE_UNKNOWN }
};

static fn_arg fn_args_project[] = {
  { "type", { .s = NULL }, 0, FN_VALTYPE_STRING },
  { NULL,   {},            0, FN_VALTYPE_UNKNOWN }
};

static fn_arg fn_args_report[] = {
  { "sumsq", { .b = 0 }, 0, FN_VALTYPE_BOOL },
  { NULL,    {},         0, FN_VALTYPE_UNKNOWN }
};

static fn_arg fn_args_resize[] = {
  { "size",  { .i = 0     }, 0, FN_VALTYPE_INT },
  { "shape", { .iv = NULL }, 0, FN_VALTYPE_INTS },
  { NULL,    {},             0, FN_VALTYPE_UNKNOWN }
};

static fn_arg fn_args_shift[] = {
  { "points", { .b = 0   }, 0, FN_VALTYPE_BOOL },
  { "sec",    { .b = 0   }, 0, FN_VALTYPE_BOOL },
  { "ppm",    { .b = 0   }, 0, FN_VALTYPE_BOOL },
  { "hz",     { .b = 0   }, 0, FN_VALTYPE_BOOL },
  { "round",  { .b = 0   }, 0, FN_VALTYPE_BOOL },
  { "amount", { .f = 0.0 }, 0, FN_VALTYPE_FLOAT },
  { NULL,     {},           0, FN_VALTYPE_UNKNOWN }
};

static fn_arg fn_args_subsamp[] = {
  { "sched", { .s = NULL }, 0, FN_VALTYPE_STRING },
  { NULL,    {},            0, FN_VALTYPE_UNKNOWN }
};

static fn_arg fn_args_symm[] = {
  { "dims", { .iv = NULL }, 0, FN_VALTYPE_INTS },
  { NULL,   {},             0, FN_VALTYPE_UNKNOWN }
};

static fn_arg fn_args_tilt[] = {
  { "angle", { .f = 0.0   }, 0, FN_VALTYPE_FLOAT },
  { "dims",  { .iv = NULL }, 0, FN_VALTYPE_INTS },
  { NULL,    {},             0, FN_VALTYPE_UNKNOWN }
};

static fn_arg fn_args_window[] = {
  { "type",   { .s = NULL }, 0, FN_VALTYPE_STRING },
  { "start",  { .f = 0.0  }, 0, FN_VALTYPE_FLOAT },
  { "end",    { .f = 1.0  }, 0, FN_VALTYPE_FLOAT },
  { "order",  { .f = 1.0  }, 0, FN_VALTYPE_FLOAT },
  { "lb",     { .f = 0.0  }, 0, FN_VALTYPE_FLOAT },
  { "invlb",  { .f = 0.0  }, 0, FN_VALTYPE_FLOAT },
  { "center", { .f = 0.0  }, 0, FN_VALTYPE_FLOAT },
  { NULL,     {},            0, FN_VALTYPE_UNKNOWN }
};

static fn_arg fn_args_zerofill[] = {
  { "times", { .i = 0 }, 0, FN_VALTYPE_INT },
  { NULL,    {},         0, FN_VALTYPE_UNKNOWN }
};

#endif /* __HXND_FN_ARGS_H__ */

