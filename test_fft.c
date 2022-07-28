#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
    float real;
    float imag;
} Complex;

int ones_32(int32_t n);
int32_t floor_log2_32(int32_t x);
int fft(Complex x[], int32_t N);
int inverse_fft(Complex x[], int32_t N);

const float sin_tb[] = {  //(PI PI/2 PI/4 PI/8 PI/16 ... PI/(2^k))
    0.000000, 1.000000, 0.707107, 0.382683, 0.195090, 0.098017, 0.049068,
    0.024541, 0.012272, 0.006136, 0.003068, 0.001534, 0.000767, 0.000383,
    0.000192, 0.000096, 0.000048, 0.000024, 0.000012, 0.000006, 0.000003};

const float cos_tb[] = {  //(PI PI/2 PI/4 PI/8 PI/16 ... PI/(2^k))
    -1.000000, 0.000000, 0.707107, 0.923880, 0.980785, 0.995185, 0.998795,
    0.999699,  0.999925, 0.999981, 0.999995, 0.999999, 1.000000, 1.000000,
    1.000000,  1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000};

#define SAMPLE_NODES (32)
#define PI (3.14159265f)

int main(void) {
    Complex x[SAMPLE_NODES];
    Complex reference[SAMPLE_NODES];

    for (int i = 0; i < SAMPLE_NODES; i++) {
        x[i].real = (0.7*sin(2*5*PI*i/SAMPLE_NODES) + sin(2*PI*3*i/SAMPLE_NODES));
        x[i].imag = 0.0f;
        // reference[i] = x[i];
    }
    for (int i = 0; i < SAMPLE_NODES; i++) {
        printf("%d\t : %.5f %.5f\n", i,x[i].real, x[i].imag);
    }
    printf("\n\n\n\n");

    fft(x, SAMPLE_NODES);

    for (int i = 0; i < SAMPLE_NODES; i++) {
        if(x[i].real<0){
            x[i].real *= -1;
        }
        if(x[i].imag<0){
            x[i].imag *= -1;
        }

        printf("%d\t : %.5f %.5f\n", i, x[i].real*2/SAMPLE_NODES, x[i].imag*2/SAMPLE_NODES);
    }
    printf("\n\n\n\n");

    inverse_fft(x, SAMPLE_NODES);

    //for (int i = 0; i < SAMPLE_NODES; i++) {
    //    printf(" %.5f %.5f\n", x[i].real, x[i].imag);
    //}

}

int ones_32(int32_t n) {
    unsigned int c = 0;
    for (c = 0; n; ++c) {
        n &= (n - 1);
    }
    return c;
}

int32_t floor_log2_32(int32_t x) {
    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);

    return (ones_32(x >> 1));
}

int fft(Complex x[], int32_t N) {
    int i, j, l, k, ip;
    static int32_t M = 0;
    static int le, le2;
    static float sR, sI, tR, tI, uR, uI;

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
    for (l = 1; l <= M; l++) { /* loop for ceil{log2(N)} */
        le = (int)(1 << l);
        le2 = (int)(le >> 1);
        uR = 1;
        uI = 0;

        k = floor_log2_32(le2);
        sR = cos_tb[k];                       // cos(PI / le2);
        sI = -sin_tb[k];                      // -sin(PI / le2)
        for (j = 1; j <= le2; j++) {          /* loop for each sub DFT */
            for (i = j - 1; i < N; i += le) { /* loop for each butterfly */
                ip = i + le2;
                tR = (x[ip].real * uR) - (x[ip].imag * uI);
                tI = (x[ip].real * uI) + (x[ip].imag * uR);
                x[ip].real = x[i].real - tR;
                x[ip].imag = x[i].imag - tI;
                x[i].real += tR;
                x[i].imag += tI;
            } /* Next i */
            tR = uR;
            uR = (tR * sR) - (uI * sI);
            uI = (tR * sI) + (uI * sR);
        } /* Next j */
    }     /* Next l */

    return 0;
}

int inverse_fft(Complex x[], int32_t N) {
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