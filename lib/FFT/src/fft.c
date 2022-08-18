#include "fft.h"

#define FLOAT_TO_INT(x) ((x) >= 0 ? (int)((x) + 0.5) : (int)((x)-0.5))

static const int sin_tb[] = {
    FLOAT_TO_INT(0.000000 * (1 << FFT_MATH_FRACTION_BITS)), //PI
    FLOAT_TO_INT(1.000000 * (1 << FFT_MATH_FRACTION_BITS)), //PI/2
    FLOAT_TO_INT(0.707107 * (1 << FFT_MATH_FRACTION_BITS)), //PI/4
    FLOAT_TO_INT(0.382683 * (1 << FFT_MATH_FRACTION_BITS)), //PI/8
    FLOAT_TO_INT(0.195090 * (1 << FFT_MATH_FRACTION_BITS)), //PI/16
    FLOAT_TO_INT(0.098017 * (1 << FFT_MATH_FRACTION_BITS)), //...
    FLOAT_TO_INT(0.049068 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(0.024541 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(0.012272 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(0.006136 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(0.003068 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(0.001534 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(0.000767 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(0.000383 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(0.000192 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(0.000096 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(0.000048 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(0.000024 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(0.000012 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(0.000006 * (1 << FFT_MATH_FRACTION_BITS)),  //...
    FLOAT_TO_INT(0.000003 * (1 << FFT_MATH_FRACTION_BITS))}; //PI/(2^K)

static const int cos_tb[] = {
    FLOAT_TO_INT(-1.000000 * (1 << FFT_MATH_FRACTION_BITS)), //PI
    FLOAT_TO_INT( 0.000000 * (1 << FFT_MATH_FRACTION_BITS)), //PI/2
    FLOAT_TO_INT( 0.707107 * (1 << FFT_MATH_FRACTION_BITS)), //PI/4
    FLOAT_TO_INT( 0.923880 * (1 << FFT_MATH_FRACTION_BITS)), //PI/8
    FLOAT_TO_INT( 0.980785 * (1 << FFT_MATH_FRACTION_BITS)), //PI/16
    FLOAT_TO_INT( 0.995185 * (1 << FFT_MATH_FRACTION_BITS)), //...
    FLOAT_TO_INT( 0.998795 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT( 0.999699 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT( 0.999925 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT( 0.999981 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT( 0.999995 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT( 0.999999 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT( 1.000000 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT( 1.000000 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT( 1.000000 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT( 1.000000 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT( 1.000000 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT( 1.000000 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT( 1.000000 * (1 << FFT_MATH_FRACTION_BITS)),
    FLOAT_TO_INT( 1.000000 * (1 << FFT_MATH_FRACTION_BITS)),  //...
    FLOAT_TO_INT( 1.000000 * (1 << FFT_MATH_FRACTION_BITS))}; //PI/(2^K)

/**
 * @brief Simple Fast Fourier Transform also known as FFT.
 *
 * @param x is a fixedpoint array of the complex datatype defined in fft.h,
 * the answer of the FFT will be in this array. The data in this
 * array will be deleted, make sure to save the data if you need it.
 *
 * @param N is the length of the array. NOTE this variable need to
 * be a number 2^k.
 *
 * @return The function returns 0, the answer is in the array.
 */
int32_t fft(Complex x[], int32_t N) {
    int i, j, l, k, ip;
    static int32_t M = 0;
    static int le, le2;
    static int sR, sI;
    static int uR, uI, tR, tI;

    M = floor_log2_32(N);
    /*
     * bit reversal sorting
     */
    l = N >> 1;
    j = l;
    ip = N - 2;
    for (i = 1; i <= ip; i++) {
        if (i < j) {
            tR = x[j].real;
            tI = x[j].imag;
            x[j].real = x[i].real;
            x[j].imag = x[i].imag;
            x[i].real = tR;
            x[i].imag = tI;
        }
        k = l;
        while (k <= j) {
            j = j - k;
            k = k >> 1;
        }
        j = j + k;
    }

    /*
     * For Loops
     */
    for (l = 1; l <= M; l++) {
        le = (int)(1 << l);
        le2 = (int)(le >> 1);
        uR = 1 << FFT_MATH_FRACTION_BITS;
        uI = 0 << FFT_MATH_FRACTION_BITS;

        k = floor_log2_32(le2);
        sR = cos_tb[k];
        sI = -sin_tb[k];
        for (j = 1; j <= le2; j++) {          /* loop for each sub DFT */
            for (i = j - 1; i < N; i += le) { /* loop for each butterfly */
                ip = i + le2;
                tR = (((long)uR * x[ip].real) >> FFT_MATH_FRACTION_BITS) -
                     (((long)uI * x[ip].imag) >> FFT_MATH_FRACTION_BITS);
                tI = (((long)uI * x[ip].real) >> FFT_MATH_FRACTION_BITS) +
                     (((long)uR * x[ip].imag) >> FFT_MATH_FRACTION_BITS);
                x[ip].real = x[i].real - tR;
                x[ip].imag = x[i].imag - tI;
                x[i].real += tR;
                x[i].imag += tI;
            } /* Next i */
            /* Calculation of twiddle factor */
            tR = uR;
            uR = (((long)tR * sR) >> FFT_MATH_FRACTION_BITS) -
                 (((long)uI * sI) >> FFT_MATH_FRACTION_BITS);
            uI = (((long)tR * sI) >> FFT_MATH_FRACTION_BITS) +
                 (((long)uI * sR) >> FFT_MATH_FRACTION_BITS);
        } /* Next j */
    }     /* Next l */

    return 0;
}

/**
 * @brief Simple inverse Fast Fourier Transform also known as IFFT. 
 * 
 * @param x is a fixedpoint array of the complex datatype defined in fft.h,
 * the answer of the IFFT will be in this array. The data in this
 * array will be deleted, make sure to save the data if you need it.
 * 
 * @param N is the length of the array. NOTE this variable need to
 * be a number 2^k.
 * 
 * @return The function returns 0, the answer is in the array. 
 */
int32_t inverse_fft(Complex x[], int32_t N) {
    int k = 0;

    for (k = 0; k <= N - 1; k++) {
        x[k].imag = -x[k].imag;
    }

    fft(x, N); /* using FFT */

    for (k = 0; k <= N - 1; k++) {
        x[k].real = x[k].real / N;
        x[k].imag = -x[k].imag / N;
    }

    return 0;
}

/**
 * @brief Returns the number of ones, bitwise in an function.
 *
 * @param n the number.
 * @return int number of ones.
 */
static int32_t ones_32(int32_t n) {
    unsigned int c = 0;
    for (c = 0; n; ++c) {
        n &= (n - 1);
    }
    return c;
}

/**
 * @brief Returns the floor(log2(x))
 *
 * @param x the number
 * @return floor log2 of x
 */
static int32_t floor_log2_32(int32_t x) {
    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);

    return (ones_32(x >> 1));
}