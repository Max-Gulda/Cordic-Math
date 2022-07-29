#pragma once

#include "stdint.h"


#define FFT_MATH_FRACTION_BITS 16

typedef struct {
    int real;
    int imag;
} Complex;

int32_t ones_32(int32_t n);
int32_t floor_log2_32(int32_t x);
int32_t fft(Complex x[], int32_t N);
int32_t inverse_fft(Complex x[], int32_t N);