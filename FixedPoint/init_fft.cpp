#include <stdio.h>

#include "fft.h"
#include "fft32.h"
#include "fft64.h"

// Maximale Länge für die Float-FFT
#define MAX_FFT 262144

// Cosinus- und Sinus-Tabelle für die 32-Bit-Fixed-Point-FFT
#define MAX_FFT32 16384
static int32_t wtable32[3 * MAX_FFT32 / 4 + 1] = {
#include "fft_table32.h"
};

// Cosinus- und Sinus-Tabelle für die 64-Bit-Fixed-Point-FFT
#define MAX_FFT64 8192
static int64_t wtable64[3 * MAX_FFT64 / 4 + 1] = {
#include "fft_table64.h"
};

void init_fft(void)
{
    printf("INIT FFT\n");
    set_twiddle_table(MAX_FFT);
    set_twiddle_table_double(MAX_FFT);

    set_twiddle_table32(wtable32, MAX_FFT32);
    set_twiddle_table64(wtable64, MAX_FFT64);
}
