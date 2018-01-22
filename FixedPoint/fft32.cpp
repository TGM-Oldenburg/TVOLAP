/*----------------------------------------------------------------------------*\
| 32-bit fixed-point radix-2 decimation-in-time fast Fourier transform (FFT)   |
|                                                                              |
| Authors:  Wiebke Lamping                                                     |
|           Simon Becker                                                       |
|           Moritz Wächtler                                                    |
|           Ben Williges                                                       |
|           Uwe Simmer                                                         |
|                                                                              |
| Date:     21.03.2011      initial version                                    |
|           26.07.2012      stride for sin/cos tables                          |
|           18.11.2012      complex and real ffts                              |
\*----------------------------------------------------------------------------*/

#include <stdio.h>
#include "fft32.h"

typedef struct {
    int nfft;
    complex32 *twiddle_factor;
    int32_t *cos_half;
} fft_table32;

static fft_table32 table32;

//-----------------------------------------------------------------------------

void set_twiddle_table32(int32_t *fft_table, int max_nfft)
{
    table32.nfft = fft_table[0];
    table32.twiddle_factor = (complex32 *) (fft_table + 1);
    table32.cos_half = fft_table + table32.nfft + 1;
}

//------------------------------------------------------------------------------

int table_get_nfft32()
{
    return table32.nfft;
}

//------------------------------------------------------------------------------

complex32 *table_get_twiddle_factor32()
{
    return table32.twiddle_factor;
}

//------------------------------------------------------------------------------

int32_t *table_get_cos_half32()
{
    return table32.cos_half;
}

//-----------------------------------------------------------------------------

void fft_core_scale32(complex32 *x, complex32 *w, int nstride, int nfft)
{
    int i, ig, j, k, l, ngroups, nbutterflies;          // 32.0
    complex32 ctemp;                                    // 1.31

    //----- Bit Reverse Section -----------------------------------------------

    j = 0;
    for (i=0; i<nfft-1; i++)
    {
        if (i<j)
        {
            ctemp.re = x[j].re;
            ctemp.im = x[j].im;

            x[j].re = x[i].re;
            x[j].im = x[i].im;

            x[i].re = ctemp.re;
            x[i].im = ctemp.im;
        }

        k = nfft >> 1;
        while (k <= j)
        {
            j -= k;
            k = k >> 1;
        }
        j += k;
    }

    //----- First Stage -------------------------------------------------------

    i = 0;
    j = 1;
    ngroups = nfft >> 1;

    for (ig=0; ig<ngroups; ig++)
    {
        ctemp.re = x[j].re >> 1;
        ctemp.im = x[j].im >> 1;

        x[j].re = (x[i].re >> 1) - ctemp.re;
        x[j].im = (x[i].im >> 1) - ctemp.im;

        x[i].re = (x[i].re >> 1) + ctemp.re;
        x[i].im = (x[i].im >> 1) + ctemp.im;

        i = i + 2;
        j = j + 2;
    }

    //----- Final log2(nfft)-1 Stages -----------------------------------------

    nbutterflies = 2;

    for (ngroups = nfft >> 2; ngroups > 0; ngroups >>= 1)
    {
        i = 0;
        j = nbutterflies;

        for (ig=0; ig<ngroups; ig++)
        {
            k = 0;
            for (l=0; l<nbutterflies; l++)
            {
                // Butterfly
                ctemp.re = (((int64_t)x[j].re * w[k].re) >> 32)
                         - (((int64_t)x[j].im * w[k].im) >> 32);
                ctemp.im = (((int64_t)x[j].im * w[k].re) >> 32)
                         + (((int64_t)x[j].re * w[k].im) >> 32);

                x[j].re = (x[i].re >> 1) - ctemp.re;
                x[j].im = (x[i].im >> 1) - ctemp.im;

                x[i].re = (x[i].re >> 1) + ctemp.re;
                x[i].im = (x[i].im >> 1) + ctemp.im;

                k += (ngroups << nstride);
                i++;
                j++;
            }

            i = i + nbutterflies;
            j = j + nbutterflies;
        }

        nbutterflies <<= 1;
    }
}

//-----------------------------------------------------------------------------

void fft_core32(complex32 *x, complex32 *w, int nstride, int nfft)
{
    int i, ig, j, k, l, ngroups, nbutterflies;          // 32.0
    complex32 ctemp;                                    // 1.31

//----- Bit Reverse Section -----------------------------------------------

    j = 0;
    for (i=0; i<nfft-1; i++)
    {
        if (i<j)
        {
            ctemp.re = x[j].re;
            ctemp.im = x[j].im;

            x[j].re = x[i].re;
            x[j].im = x[i].im;

            x[i].re = ctemp.re;
            x[i].im = ctemp.im;
        }

        k = nfft >> 1;

        while (k <= j)
        {
            j -= k;
            k = k >> 1;
        }
        j += k;
    }

    //----- First Stage -------------------------------------------------------

    i = 0;
    j = 1;
    ngroups = nfft >> 1;

    for (ig=0; ig<ngroups; ig++)
    {
        ctemp.re = x[j].re;
        ctemp.im = x[j].im;

        x[j].re = x[i].re - ctemp.re;
        x[j].im = x[i].im - ctemp.im;

        x[i].re = x[i].re + ctemp.re;
        x[i].im = x[i].im + ctemp.im;

        i = i + 2;
        j = j + 2;
    }

    //----- Final log2(nfft)-1 Stages -----------------------------------------

    nbutterflies = 2;

    for (ngroups = nfft >> 2; ngroups > 0; ngroups >>= 1)
    {
        i = 0;
        j = nbutterflies;

        for (ig=0; ig<ngroups; ig++)
        {
            k = 0;
            for (l=0; l<nbutterflies; l++)
            {
                // Butterfly
                ctemp.re = ((((int64_t)x[j].re * w[k].re) >> 32)
                          - (((int64_t)x[j].im * w[k].im) >> 32)) << 1;
                ctemp.im = ((((int64_t)x[j].im * w[k].re) >> 32)
                          + (((int64_t)x[j].re * w[k].im) >> 32)) << 1;

                x[j].re = x[i].re - ctemp.re;
                x[j].im = x[i].im - ctemp.im;

                x[i].re = x[i].re + ctemp.re;
                x[i].im = x[i].im + ctemp.im;

                k += (ngroups << nstride);
                i++;
                j++;
            }

            i = i + nbutterflies;
            j = j + nbutterflies;
        }
        nbutterflies <<= 1;
    }
}

//-----------------------------------------------------------------------------

void rfft32(int32_t *input, complex32 *spectrum, int scale, int n)
{
    int i, j, k, nfft, nstride;                         // 32.0
    int32_t tr, ti, rs, is, rd, id, rp, ip, ci, cj;     // 1.31
    complex32 *x, *w;                                   // 1.31
    int32_t *cos2table;                                 // 1.31

    nfft = n >> 1;

    if (nfft == 0 || input == NULL || spectrum == NULL)
        return;

    // Table stride
    k = table32.nfft;
    nstride = -1;
    i = 0;
    while (k)
    {
        if (k == nfft)
        {
            nstride = i;
            break;
        }
        k = k >> 1;
        i++;
    }

    if (nstride < 0)
        return;

    // Copy memory if not in-place
    if (input != (int32_t *) spectrum)
    {
        for (i=0; i<n/2; i++)
        {
            spectrum[i].re = input[2*i + 0];
            spectrum[i].im = input[2*i + 1];
        }
    }

    x = spectrum;
    w = table32.twiddle_factor;
    cos2table = table32.cos_half;

    if (scale)
        fft_core_scale32(x, w, nstride, nfft);
    else
        fft_core32(x, w, nstride, nfft);

    //----- Half Length Postprocessing ----------------------------------------

    tr = x[0].re;
    ti = x[0].im;

    x[0].re = (tr >> 1) + (ti >> 1);
    x[0].im = 0;

    x[nfft].re = (tr >> 1) - (ti >> 1);
    x[nfft].im = 0;

    x[nfft/2].re = x[nfft/2].re >> 1;
    x[nfft/2].im = -x[nfft/2].im >> 1;

    for (i=1; i<nfft/2; i++)
    {
        j = nfft - i;

        rs = (x[i].re >> 1) + (x[j].re >> 1);
        rd = (x[j].re >> 1) - (x[i].re >> 1);
        is = (x[i].im >> 1) + (x[j].im >> 1);
        id = (x[i].im >> 1) - (x[j].im >> 1);

        ci = cos2table[i << nstride];
        cj = cos2table[(nfft/2-i) << nstride];

        rp = (((int64_t) is * ci) >> 31)
           + (((int64_t) rd * cj) >> 31);
        ip = (((int64_t) rd * ci) >> 31)
           - (((int64_t) is * cj) >> 31);

        x[i].re = (rp >> 1) + (rs >> 1);
        x[j].re = (rs >> 1) - (rp >> 1);

        x[i].im = (ip >> 1) + (id >> 1);
        x[j].im = (ip >> 1) - (id >> 1);
    }

    if (!scale)
    {
        for (i=0; i<nfft; i++)
        {
            x[i].re <<= 1;
            x[i].im <<= 1;
        }
    }
}


//-----------------------------------------------------------------------------

void irfft32(complex32 *spectrum, int32_t *output, int scale, int n)
{
    int i, j, k, nfft, nstride;                         // 32.0
    int32_t t0, tn, rs, is, rd, id, rp, ip, ci, cj;     // 1.31
    complex32 *x, *w;                                   // 1.31
    int32_t *cos2table;                                 // 1.31

    nfft = n >> 1;

    if (nfft == 0 || spectrum == NULL || output == NULL)
        return;

    // Table stride
    k = table32.nfft;
    nstride = -1;
    i = 0;
    while (k)
    {
        if (k == nfft)
        {
            nstride = i;
            break;
        }
        k = k >> 1;
        i++;
    }

    if (nstride < 0)
        return;

    x = (complex32 *) output;
    w = table32.twiddle_factor;
    cos2table = table32.cos_half;

    //----- Half Length Preprocessing -----------------------------------------

    t0 = spectrum[0].re >> 1;
    tn = spectrum[nfft].re >> 1;

    x[0].re = (t0 + tn);
    x[0].im = (t0 - tn);

    x[nfft/2].re = spectrum[nfft/2].re;
    x[nfft/2].im = -spectrum[nfft/2].im;

    for (i=1; i<(nfft>>1); i++)
    {
        j = nfft - i;

        rs = (spectrum[i].re + spectrum[j].re) >> 1;
        rd = (spectrum[i].re - spectrum[j].re) >> 1;
        is = (spectrum[i].im + spectrum[j].im) >> 1;
        id = (spectrum[i].im - spectrum[j].im) >> 1;

        ci = cos2table[i << nstride];
        cj = cos2table[(nfft/2-i) << nstride];

        rp = (((int64_t) is * ci) >> 31)
           + (((int64_t) rd * cj) >> 31);
        ip = (((int64_t) rd * ci) >> 31)
           - (((int64_t) is * cj) >> 31);

        x[i].re = rp + rs;
        x[j].re = rs - rp;

        x[i].im = ip - id;
        x[j].im = ip + id;
    }

    if (scale)
        fft_core_scale32(x, w, nstride, nfft);
    else
        fft_core32(x, w, nstride, nfft);

    if (!scale)
    {
        for (i=0; i<nfft; i++)
        {
            x[i].re <<= 1;
            x[i].im <<= 1;
        }
    }
}

//-----------------------------------------------------------------------------

void cfft32(complex32 *x, int scale, int nfft)
{
    int i, k, nstride;

    if (nfft == 0 || x == NULL)
        return;

    // Table stride
    k = table32.nfft;
    nstride = -1;
    i = 0;
    while (k)
    {
        if (k == nfft)
        {
            nstride = i;
            break;
        }
        k = k >> 1;
        i++;
    }

    if (nstride < 0)
        return;

    if (scale)
        fft_core_scale32(x, table32.twiddle_factor, nstride, nfft);
    else
        fft_core32(x, table32.twiddle_factor, nstride, nfft);
}

//-----------------------------------------------------------------------------

void icfft32(complex32 *x, int scale, int nfft)
{
    int i, k, nstride;

    if (nfft == 0 || x == NULL)
        return;

    // Table stride
    k = table32.nfft;
    nstride = -1;
    i = 0;
    while (k)
    {
        if (k == nfft)
        {
            nstride = i;
            break;
        }
        k = k >> 1;
        i++;
    }

    if (nstride < 0)
        return;

    for (i=0; i<nfft; i++)
        x[i].im = -x[i].im;

    if (scale)
        fft_core_scale32(x, table32.twiddle_factor, nstride, nfft);
    else
        fft_core32(x, table32.twiddle_factor, nstride, nfft);

    for (i=0; i <nfft; i++)
        x[i].im = -x[i].im;
}

//--------------------- Licence -----------------------------------------------
// Copyright (c) <2011-2012> Wiebke Lamping, Simon Becker, Moritz Wächtler, 
// Ben Williges, Uwe Simmer, Institute for Hearing Technology and Audiology 
// Jade University of Applied Sciences Oldenburg 
// Permission is hereby granted, free of charge, to any person obtaining 
// a copy of this software and associated documentation files 
// (the "Software"), to deal in the Software without restriction, including 
// without limitation the rights to use, copy, modify, merge, publish, 
// distribute, sublicense, and/or sell copies of the Software, and to 
// permit persons to whom the Software is furnished to do so, subject 
// to the following conditions:
// The above copyright notice and this permission notice shall be included 
// in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
