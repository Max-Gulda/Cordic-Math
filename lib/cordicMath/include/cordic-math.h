/* Copyright (c) 2022 Max Gulda, KTH

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE. */
#include "gd32vf103.h"


#define PI 804
#define EULER 696
#define CORDIC_GAIN 156
#define DECIMAL_TO_FP 256


int LUT_CORDIC_ATAN[15] = { 11520,  /* 45.000    degrees */
                            6801,   /* 26.566    degrees */
                            3593,   /* 26.566    degrees */
                            1824,   /* 14.035    degrees */
                            916,    /* 3.578     degrees */
                            458,    /* 1.789     degrees */
                            229,    /* 0.894     degrees */
                            115,    /* 0.449     degrees */
                            57,     /* 0.223     degrees */
                            28,     /* 0.109     degrees */
                            14,     /* 0.055     degrees */
                            7,      /* 0.027     degrees */
                            4,      /* 0.016     degrees */
                            2,      /* 0.008     degrees */
                            1};     /* 0.004     degrees */
 
int LUT_CORDIC_ATANH[14] = {8057,   /* 31.473    degrees */
                            3746,   /* 14.633    degrees */
                            1843,   /* 7.199     degrees */
                            918,    /* 3.586     degrees */
                            459,    /* 1.793     degrees */
                            229,    /* 0.895     degrees */
                            115,    /* 0.449     degrees */
                            52,     /* 0.203     degrees */
                            29,     /* 0.113     degrees */
                            14,     /* 0.055     degrees */
                            7,      /* 0.027     degrees */
                            4,      /* 0.016     degrees */
                            2,      /* 0.008     degrees */
                            1};     /* 0.004     degrees */

int32_t cordic_atan(int32_t y, int32_t x);
int32_t cordic_hypotenuse(int32_t y, int32_t x);
int32_t cordic_cos(int32_t theta);
int32_t cordic_sin(int32_t theta);
int32_t cordic_asin(int32_t yInput);
int32_t cordic_acos(int32_t xInput);
int32_t cordic_tan(int32_t degree);
int32_t cordic_sqrt(int32_t x);
int32_t abs(int32_t input);
int32_t isEven(int32_t input);
int32_t isOdd(int32_t input);
int32_t to_degree(int32_t input);
int32_t to_radians(int32_t input);
int32_t cordic_arctanh(int32_t y, int32_t x);
int32_t cordic_ln(int32_t input);
int32_t cordic_arccosh(int32_t x);
int32_t cordic_arcsinh(int32_t y);
int32_t cordic_sinh(int32_t theta);
int32_t cordic_cosh(int32_t theta);
int32_t cordic_tanh(int32_t theta);
int32_t cordic_exp(int32_t exponent);
int32_t cordic_pow(int32_t base, int32_t exponent);