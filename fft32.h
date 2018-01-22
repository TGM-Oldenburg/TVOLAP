#ifndef _FFT32
#define _FFT32

/*
32-bit fixed-point radix-2 decimation-in-time FFT

Example:

#include "fft32.h"
#define NFFT 512

#define MAX_FFT32 16384
int32_t wtable32[3*MAX_FFT32/4+1] = {
#include "FFTtable32.h"
};

complex32 spectrum32[NFFT/2+1];
int32_t input32[NFFT];
int32_t output32[NFFT];

// initialization (once)
set_twiddle_table32(wtable32, MAX_FFT32);

// fft computation
fft32(input32, spectrum32, 1, NFFT);

// ifft computation
ifft32(spectrum32, output32, 0, NFFT);
*/

#include "complex32.h"

void set_twiddle_table32(int32_t *fft_table, int max_nfft);
void rfft32(int32_t *input, complex32 *spectrum, int scale, int nfft);
void irfft32(complex32 *spectrum, int32_t *output, int scale, int nfft);
void cfft32(complex32 *x, int scale, int nfft);
void icfft32(complex32 *x, int scale, int nfft);

int table_get_nfft32();
complex32 *table_get_twiddle_factor32();
int32_t *table_get_cos_half32();

#endif
