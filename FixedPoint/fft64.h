#ifndef _FFT_64
#define _FFT_64

/*
64-bit fixed-point radix-2 decimation-in-time FFT

Example:

#include "fft64.h"
#define NFFT 512

#define MAX_FFT64 16384
int64_t wtable64[3*MAX_FFT64/4+1] = {
#include "FFTtable64.h"
};

complex64 spectrum64[NFFT/2+1];
int64_t input64[NFFT];
int64_t output64[NFFT];

// initialization (once)
set_twiddle_table64(wtable64, MAX_FFT64);

// fft computation
rfft64(input64, spectrum64, 1, NFFT);

// ifft computation
irfft64(spectrum64, output64, 0, NFFT);
*/

#include "complex64.h"

void set_twiddle_table64(int64_t *fft_table, int max_nfft);
void rfft64(int64_t *input, complex64 *spectrum, int scale, int nfft);
void irfft64(complex64 *spectrum, int64_t *output, int scale, int nfft);
void cfft64(complex64 *x, int scale, int nfft);
void icfft64(complex64 *x, int scale, int nfft);

int table_get_nfft64();
complex64 *table_get_twiddle_factor64();
int64_t *table_get_cos_half64();

#endif
