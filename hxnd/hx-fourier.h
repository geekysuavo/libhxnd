
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
#ifndef __HXND_HX_FOURIER_H__
#define __HXND_HX_FOURIER_H__

/* define constants for forward and reverse fft operations.
 */
#define HX_FFT_FORWARD   1.0
#define HX_FFT_REVERSE  -1.0

/* function declarations: */

int hx_ispow2 (unsigned int value);

unsigned int hx_prevpow2 (unsigned int value);

unsigned int hx_nextpow2 (unsigned int value);

int hx_array_fft1d (hx_array *y, int d, real dir,
                    hx_scalar *w, hx_scalar *swp);

int hx_array_fftfn (hx_array *x, int d, int k, real dir);

#define hx_array_fft(x, d, k) \
  hx_array_fftfn(x, d, k, HX_FFT_FORWARD)

#define hx_array_ifft(x, d, k) \
  hx_array_fftfn(x, d, k, HX_FFT_REVERSE)

int hx_array_ht (hx_array *x, int d, int k);

int hx_array_fshift (hx_array *x, int d, int k, real amount);

#endif /* __HXND_HX_FOURIER_H__ */

