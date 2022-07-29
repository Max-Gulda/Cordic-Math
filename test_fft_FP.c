#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define FFT_MATH_FRACTION_BITS 8
#define SAMPLE_NODES (16)
#define PI (3.14159265f)

#define FLOAT_TO_INT(x) ((x) >= 0 ? (int)((x) + 0.5) : (int)((x)-0.5))

typedef struct {
    int real;
    int imag;
} Complex;

int32_t ones_32(int32_t n);
int32_t floor_log2_32(int32_t x);
int32_t fft(Complex x[], int32_t N);
int32_t inverse_fft(Complex x[], int32_t N);

static const int sin_tb[] = {  //(PI PI/2 PI/4 PI/8 PI/16 ... PI/(2^k))
                        FLOAT_TO_INT( 0.000000 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 1.000000 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 0.707107 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 0.382683 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 0.195090 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 0.098017 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 0.049068 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 0.024541 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 0.012272 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 0.006136 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 0.003068 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 0.001534 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 0.000767 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 0.000383 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 0.000192 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 0.000096 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 0.000048 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 0.000024 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 0.000012 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 0.000006 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 0.000003 * (1 << FFT_MATH_FRACTION_BITS))};

static const int cos_tb[] = {  //(PI PI/2 PI/4 PI/8 PI/16 ... PI/(2^k))
                        FLOAT_TO_INT(-1.000000 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 0.000000 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 0.707107 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 0.923880 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 0.980785 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 0.995185 * (1 << FFT_MATH_FRACTION_BITS)),
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
                        FLOAT_TO_INT( 1.000000 * (1 << FFT_MATH_FRACTION_BITS)),
                        FLOAT_TO_INT( 1.000000 * (1 << FFT_MATH_FRACTION_BITS))};



int main(void) {
    Complex x[SAMPLE_NODES];
    Complex reference[SAMPLE_NODES];

    for (int i = 0; i < SAMPLE_NODES; i++) {
        x[i].real = (0.7*sin(2*5*PI*i/SAMPLE_NODES) + sin(2*PI*3*i/SAMPLE_NODES))* (1 << FFT_MATH_FRACTION_BITS);
        x[i].imag = 0;
        // reference[i] = x[i];
    }
    for (int i = 0; i < SAMPLE_NODES; i++) {
        printf("%d\t : %.5f %.5f\n", i,x[i].real/pow(2,FFT_MATH_FRACTION_BITS), x[i].imag/pow(2,FFT_MATH_FRACTION_BITS));
    }
    printf("\n\n\n\n");

    fft(x, SAMPLE_NODES);

    for (int i = 0; i < SAMPLE_NODES; i++) {
        printf("%d\t : %.5f %.5f\n", i, (x[i].real*2/SAMPLE_NODES)/pow(2,FFT_MATH_FRACTION_BITS), (x[i].imag*2/SAMPLE_NODES)/pow(2,FFT_MATH_FRACTION_BITS));
    }
    printf("\n\n\n\n");
    

    inverse_fft(x, SAMPLE_NODES);

    for (int i = 0; i < SAMPLE_NODES; i++) {
        printf("%d\t : %.5f %.5f\n", i , x[i].real/pow(2,FFT_MATH_FRACTION_BITS), x[i].imag/pow(2,FFT_MATH_FRACTION_BITS));
    }

}

int fft(Complex x[], int32_t N)
{
	int i,j,l,k,ip;
	static int32_t M = 0;
	static int le,le2;
	static int sR,sI,tR,tI,uR,uI;

	M = floor_log2_32(N);
	/*
	 * bit reversal sorting
	 */
	l = N >> 1;
	j = l;
    ip = N-2;
    for (i=1; i<=ip; i++) {
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
	for (l=1; l<=M; l++) {   /* loop for ceil{log2(N)} */
		//le = (int)pow(2,l);
		le  = (int)(1 << l) ;
		le2 = (int)(le >> 1) ;
		uR = 1 << FFT_MATH_FRACTION_BITS;
		uI = 0 << FFT_MATH_FRACTION_BITS;

        k = floor_log2_32(le2);
        sR = cos_tb[k]; //cos(PI / le2);
        sI = -sin_tb[k];  // -sin(PI / le2)
		for (j=1; j<=le2; j++) {   /* loop for each sub DFT */
			//jm1 = j - 1;
			for (i=j-1; i<N; i+=le) {  /* loop for each butterfly */
				ip = i + le2;
				tR = ((long)(x[ip].real * uR) >> FFT_MATH_FRACTION_BITS) - (((long)x[ip].imag * uI) >> FFT_MATH_FRACTION_BITS);
				tI = (((long)x[ip].real * uI) >> FFT_MATH_FRACTION_BITS) + (((long)x[ip].imag * uR) >> FFT_MATH_FRACTION_BITS);
				x[ip].real = x[i].real - tR;
				x[ip].imag = x[i].imag - tI;
				x[i].real += tR;
				x[i].imag += tI;
			}  /* Next i */
			tR = uR;
			uR = (((long)tR * sR) >> FFT_MATH_FRACTION_BITS) - (((long)uI * sI) >> FFT_MATH_FRACTION_BITS);
			uI = (((long)tR * sI) >> FFT_MATH_FRACTION_BITS) + (((long)uI * sR) >> FFT_MATH_FRACTION_BITS);
		} /* Next j */
	} /* Next l */

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