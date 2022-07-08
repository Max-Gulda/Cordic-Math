#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define CORDIC_MATH_FRACTION_BITS 16

#define FLOAT_TO_INT(x) ((x) >= 0 ? (int)((x) + 0.5) : (int)((x)-0.5))

static const uint32_t EULER = FLOAT_TO_INT(2.71828182846 * (1 << CORDIC_MATH_FRACTION_BITS));
static const uint32_t CORDIC_GAIN = FLOAT_TO_INT(0.607253 * (1 << CORDIC_MATH_FRACTION_BITS));
static const uint32_t CORDIC_GAIN_HYPERBOLIC_VECTOR = FLOAT_TO_INT(0.82816 * (1 << CORDIC_MATH_FRACTION_BITS));
static const uint32_t CORDIC_GAIN_HYPERBOLIC_CIRCULAR = FLOAT_TO_INT(1.64676 * (1 << CORDIC_MATH_FRACTION_BITS));
static const uint32_t DECIMAL_TO_FP = (1 << CORDIC_MATH_FRACTION_BITS);
static const uint32_t PI = FLOAT_TO_INT(3.14159265359 * (1 << CORDIC_MATH_FRACTION_BITS));
static const uint32_t ONE_EIGHTY_DIV_PI = FLOAT_TO_INT((180 / 3.14159265359) * (1 << CORDIC_MATH_FRACTION_BITS));
static const uint32_t ONE_DIV_CORDIC_GAIN_HYPERBOLIC = FLOAT_TO_INT((1.0 / 0.82816) * (1 << CORDIC_MATH_FRACTION_BITS));


static const uint32_t LUT_CORDIC_ATAN[15] =  {FLOAT_TO_INT(45.0000 * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 45.000    degrees */
                                              FLOAT_TO_INT(26.5651 * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 26.566    degrees */
                                              FLOAT_TO_INT(14.0362 * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 26.566    degrees */
                                              FLOAT_TO_INT(7.1250  * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 14.035    degrees */
                                              FLOAT_TO_INT(3.5763  * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 3.578     degrees */
                                              FLOAT_TO_INT(1.7899  * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 1.789     degrees */
                                              FLOAT_TO_INT(0.8952  * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 0.894     degrees */
                                              FLOAT_TO_INT(0.4476  * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 0.449     degrees */
                                              FLOAT_TO_INT(0.2238  * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 0.223     degrees */
                                              FLOAT_TO_INT(0.1119  * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 0.109     degrees */
                                              FLOAT_TO_INT(0.0560  * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 0.055     degrees */
                                              FLOAT_TO_INT(0.0280  * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 0.027     degrees */
                                              FLOAT_TO_INT(0.0140  * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 0.016     degrees */
                                              FLOAT_TO_INT(0.0070  * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 0.008     degrees */
                                              FLOAT_TO_INT(0.0035  * (1 << CORDIC_MATH_FRACTION_BITS))}; /* 0.004     degrees */

static const uint32_t LUT_CORDIC_ATANH[14] = {FLOAT_TO_INT(31.4729 * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 31.473    degrees */
                                              FLOAT_TO_INT(14.6341 * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 14.633    degrees */
                                              FLOAT_TO_INT(7.1996  * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 7.199     degrees */
                                              FLOAT_TO_INT(3.5857  * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 3.586     degrees */
                                              FLOAT_TO_INT(1.7911  * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 1.793     degrees */
                                              FLOAT_TO_INT(0.8953  * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 0.895     degrees */
                                              FLOAT_TO_INT(0.4476  * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 0.449     degrees */
                                              FLOAT_TO_INT(0.2238  * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 0.203     degrees */
                                              FLOAT_TO_INT(0.1119  * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 0.113     degrees */
                                              FLOAT_TO_INT(0.0560  * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 0.055     degrees */
                                              FLOAT_TO_INT(0.0280  * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 0.027     degrees */
                                              FLOAT_TO_INT(0.0140  * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 0.016     degrees */
                                              FLOAT_TO_INT(0.0070  * (1 << CORDIC_MATH_FRACTION_BITS)),  /* 0.008     degrees */
                                              FLOAT_TO_INT(0.0035  * (1 << CORDIC_MATH_FRACTION_BITS))}; /* 0.004     degrees */

typedef struct {
	int x;
	int y;
	int theta;
    int r;

} Coordinates;

int32_t cordic_atan(int32_t y, int32_t x);
int32_t cordic_hypotenuse(int32_t y, int32_t x);
int32_t cordic_cos(int32_t theta);
int32_t cordic_sin(int32_t theta);
int32_t cordic_asin(int32_t yInput);
int32_t cordic_acos(int32_t xInput);
int32_t cordic_tan(int32_t theta);
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
int32_t cordic_rectangular_polar(Coordinates *input);
int32_t cordic_polar_rectangular(Coordinates *input);

int main(void) {
    /*Arctan*/
    printf("\nArctan  : %f Cordic-Arctan : %f\n\n",
           atan2(2.3, 4.2) * (180 / M_PI),
           cordic_atan((1 << CORDIC_MATH_FRACTION_BITS) * 2.3,
                       (1 << CORDIC_MATH_FRACTION_BITS) * 4.2) /
               pow(2, CORDIC_MATH_FRACTION_BITS));

    /*Hypotenuse*/
    printf("Hyp     : %f Cordic-Hyp     : %f\n\n",
           sqrt((2.3 * 2.3) + (4.2 * 4.2)),
           cordic_hypotenuse((1 << CORDIC_MATH_FRACTION_BITS) * 2.3,
                             (1 << CORDIC_MATH_FRACTION_BITS) * 4.2) /
               pow(2, CORDIC_MATH_FRACTION_BITS));

    /*Cosinus*/
    printf("Cosinus : %f Cordic-Cosinus : %f\n\n", cos(372*0.0174532925),
           cordic_cos(372 * (1 << CORDIC_MATH_FRACTION_BITS)) /
               pow(2, CORDIC_MATH_FRACTION_BITS));

    /*Sinus*/
    printf("Sinus   : %f  Cordic-Sinus  : %f\n\n", sin(372*0.0174532925), //90 - 270
           cordic_sin(372 * (1 << CORDIC_MATH_FRACTION_BITS)) /
               pow(2, CORDIC_MATH_FRACTION_BITS));

    /*Arc Sin*/
    printf("Arcsin  : %f Cordic-Arcsin : %f\n\n", asin(0.43) * (180 / M_PI),
           cordic_asin(0.43 * (1 << CORDIC_MATH_FRACTION_BITS)) /
               pow(2, CORDIC_MATH_FRACTION_BITS));

    /*Arc Cos*/
    printf("Arccos  : %f Cordic-Arccos : %f\n\n", acos(0.43) * (180 / M_PI),
           cordic_acos(0.43 * (1 << CORDIC_MATH_FRACTION_BITS)) /
               pow(2, CORDIC_MATH_FRACTION_BITS));

    /*Tangens*/
    printf("Tangen  : %f Cordic-Tangen  : %f\n\n", tan(25 * 0.0174532925),
           cordic_tan(25 * (1 << CORDIC_MATH_FRACTION_BITS)) /
               pow(2, CORDIC_MATH_FRACTION_BITS));

    /*Square root*/
    printf("Sqrt    : %f Cordic-Sqrt    : %f\n\n", sqrt(89.3),
           cordic_sqrt(89.3 * (1 << CORDIC_MATH_FRACTION_BITS)) /
               pow(2, CORDIC_MATH_FRACTION_BITS));

    /*Arctan H*/
    printf("Atanh   : %f Cordic-Atanh  : %f\n\n",
           atanh(0.2411566291) * (180 / M_PI),
           cordic_arctanh((1 << CORDIC_MATH_FRACTION_BITS) * 0.635590,
                          (1 << CORDIC_MATH_FRACTION_BITS) * 2.635590) /
               pow(2, CORDIC_MATH_FRACTION_BITS));

    /*Ln*/
    printf("Ln      : %f Cordic-Ln      : %f\n\n", log(89.3),
           cordic_ln(89.3 * (1 << CORDIC_MATH_FRACTION_BITS)) /
               pow(2, CORDIC_MATH_FRACTION_BITS));

    /*Arccosh*/
    printf("Arccosh : %f Cordic-Acosh  : %f\n\n", acosh(2.3) * (180 / M_PI),
           cordic_arccosh(2.3 * (1 << CORDIC_MATH_FRACTION_BITS)) /
               pow(2, CORDIC_MATH_FRACTION_BITS));

    /*Arcsinh*/
    printf("Arcsinh : %f Cordic-Asinh  : %f\n\n", asinh(2.3) * (180 / M_PI),
           cordic_arcsinh(2.3 * (1 << CORDIC_MATH_FRACTION_BITS)) /
               pow(2, CORDIC_MATH_FRACTION_BITS));

    /*Sinh*/
    printf("Sinush  : %f  Cordic-Sinush : %f\n\n", sinh(2.3 * M_PI / 180),
           cordic_sinh(2.3 * (1 << CORDIC_MATH_FRACTION_BITS)) /
               pow(2, CORDIC_MATH_FRACTION_BITS));

    /*Cosh*/
    printf("Coshh   : %f  Cordic-Cosh   : %f\n\n", cosh(24.7 * M_PI / 180),
           cordic_cosh(24.7 * (1 << CORDIC_MATH_FRACTION_BITS)) /
               pow(2, CORDIC_MATH_FRACTION_BITS));

    /*Tanh*/
    printf("Tanh    : %f  Cordic-Tanh   : %f\n\n", tanh(23.2 * M_PI / 180),
           cordic_tanh(23.2 * (1 << CORDIC_MATH_FRACTION_BITS)) /
               pow(2, CORDIC_MATH_FRACTION_BITS));

    /*Exp*/
    printf("Exp     : %f Cordic-Exp    : %f\n\n", exp(3.39),
           cordic_exp(3.39 * (1 << CORDIC_MATH_FRACTION_BITS)) /
               pow(2, CORDIC_MATH_FRACTION_BITS));

    /*Pow*/
    printf("Pow     : %f Cordic-Pow     : %f\n\n", log(89.3),
           cordic_ln(89.3 * (1 << CORDIC_MATH_FRACTION_BITS)) /
               pow(2, CORDIC_MATH_FRACTION_BITS));
    
    /*Rectangular to polar*/
    Coordinates cor;
    cor.x = 1 * (1 << CORDIC_MATH_FRACTION_BITS);
    cor.y = 1 * (1 << CORDIC_MATH_FRACTION_BITS);
    cordic_rectangular_polar(&cor);
    printf("Rectangular to polar\n r = %f   theta = %f\n",cor.r/pow(2, CORDIC_MATH_FRACTION_BITS),cor.theta/pow(2, CORDIC_MATH_FRACTION_BITS));
    
    /*Polar to rectangular*/
    cordic_polar_rectangular(&cor);
    printf("Polar to rectangular\n x = %f   y = %f\n",cor.x/pow(2, CORDIC_MATH_FRACTION_BITS),cor.y/pow(2, CORDIC_MATH_FRACTION_BITS));

}

int32_t cordic_atan(int32_t y, int32_t x) {
    int sumAngle = 0, tempX;
    if (x < 0) {
        x = -x;
        y = -y;
    }
    for (int i = 0; i < 15; i++) {
        tempX = x;
        if (y > 0) {
            /* Rotate clockwise */
            x += (y >> i);
            y -= (tempX >> i);
            sumAngle += LUT_CORDIC_ATAN[i];
        } else {
            /* Rotate counterclockwise */
            x -= (y >> i);
            y += (tempX >> i);
            sumAngle -= LUT_CORDIC_ATAN[i];
        }
    }
    return sumAngle;
}

int32_t cordic_hypotenuse(int32_t y, int32_t x) {
    int tempX;
    x = abs(x);
    y = abs(y);

    for (int i = 0; i < 15; i++) {
        tempX = x;
        if (y > 0) {
            /* Rotate clockwise */
            x += (y >> i);
            y -= (tempX >> i);
        } else {
            /* Rotate counterclockwise */
            x -= (y >> i);
            y += (tempX >> i);
        }
    }

    return ((long)x  * CORDIC_GAIN) >> CORDIC_MATH_FRACTION_BITS ;
}

int32_t cordic_cos(int32_t theta) {
    int x = CORDIC_GAIN, y = 0, sumAngle = 0, tempX;

    theta %= (360 << CORDIC_MATH_FRACTION_BITS);
    
    if (theta > (90 << CORDIC_MATH_FRACTION_BITS)) {
        sumAngle = 180 << CORDIC_MATH_FRACTION_BITS;
    }
    if (theta > (270 << CORDIC_MATH_FRACTION_BITS)) {
        sumAngle = 360 << CORDIC_MATH_FRACTION_BITS;
    }

    for (int i = 0; i < 15; i++) {
        tempX = x;
        if (theta > sumAngle) {
            /* Rotate counter clockwise */
            x -= (y >> i);
            y += (tempX >> i);
            sumAngle += LUT_CORDIC_ATAN[i];
        } else {
            /* Rotate clockwise */
            x += (y >> i);
            y -= (tempX >> i);
            sumAngle -= LUT_CORDIC_ATAN[i];
        }
    }
    if (theta > (90 << CORDIC_MATH_FRACTION_BITS) &&
        theta < (270 << CORDIC_MATH_FRACTION_BITS)) {

        x = -x;
    }
    return x;
}

int32_t cordic_sin(int32_t theta) {
    int x = CORDIC_GAIN, y = 0, sumAngle = 0, tempX;

    theta %= (360 << CORDIC_MATH_FRACTION_BITS);

    if (theta > (90 << CORDIC_MATH_FRACTION_BITS)) {
        sumAngle = 180 << CORDIC_MATH_FRACTION_BITS;
    }
    if (theta > (270 << CORDIC_MATH_FRACTION_BITS)) {
        sumAngle = 360 << CORDIC_MATH_FRACTION_BITS;
    }

    for (int i = 0; i < 15; i++) {
        tempX = x;
        if (theta > sumAngle) {
            /* Rotate counter clockwise */
            x -= (y >> i);
            y += (tempX >> i);
            sumAngle += LUT_CORDIC_ATAN[i];
        } else {
            /* Rotate clockwise */
            x += (y >> i);
            y -= (tempX >> i);
            sumAngle -= LUT_CORDIC_ATAN[i];
        }

    }

    if (theta > (90 << CORDIC_MATH_FRACTION_BITS) &&
        theta < (270 << CORDIC_MATH_FRACTION_BITS)) {

        y = -y;
    }

    return y;
}

int32_t cordic_asin(int32_t input) {
    int x = CORDIC_GAIN, y = 0, sumAngle = 0, tempX,
        ninety = (90 << CORDIC_MATH_FRACTION_BITS);

    for (int i = 0; i < 15; i++) {
        tempX = x;
        if (y < input) {
            /* Rotate counter clockwise */
            x -= (y >> i);
            y += (tempX >> i);
            sumAngle += LUT_CORDIC_ATAN[i];
        } else {
            /* Rotate clockwise */
            x += (y >> i);
            y -= (tempX >> i);
            sumAngle -= LUT_CORDIC_ATAN[i];
        }
    }
    if (sumAngle < -ninety) {
        sumAngle = -ninety;
    } else if (sumAngle > ninety) {
        sumAngle = ninety;
    }
    return sumAngle;
}

int32_t cordic_acos(int32_t xInput) {
    int x = 0, y = CORDIC_GAIN, sumAngle = 90 << CORDIC_MATH_FRACTION_BITS,
        tempX;

    for (int i = 0; i < 15; i++) {
        tempX = x;
        if (x > xInput) {
            /* Rotate counter clockwise */
            x -= (y >> i);
            y += (tempX >> i);
            sumAngle += LUT_CORDIC_ATAN[i];
        } else {
            /* Rotate clockwise */
            x += (y >> i);
            y -= (tempX >> i);
            sumAngle -= LUT_CORDIC_ATAN[i];
        }
    }
    if (sumAngle > 180 * DECIMAL_TO_FP) {
        sumAngle = 180 * DECIMAL_TO_FP;
    } else if (sumAngle < 0) {
        sumAngle = 0;
    }
    return sumAngle;
}

//int32_t cordic_tan(int32_t theta) {
//    return ((cordic_sin(theta) << CORDIC_MATH_FRACTION_BITS) / cordic_cos(theta));
//}

int32_t cordic_tan(int32_t theta){
    int x = CORDIC_GAIN, y = 0, sumAngle = 0, tempX;

    theta %= (360 << CORDIC_MATH_FRACTION_BITS);

    if (theta > (90 << CORDIC_MATH_FRACTION_BITS)) {
        sumAngle = 180 << CORDIC_MATH_FRACTION_BITS;
    }
    if (theta > (270 << CORDIC_MATH_FRACTION_BITS)) {
        sumAngle = 360 << CORDIC_MATH_FRACTION_BITS;
    }

    for (int i = 0; i < 15; i++) {
        tempX = x;
        if (theta > sumAngle) {
            /* Rotate counter clockwise */
            x -= (y >> i);
            y += (tempX >> i);
            sumAngle += LUT_CORDIC_ATAN[i];
        } else {
            /* Rotate clockwise */
            x += (y >> i);
            y -= (tempX >> i);
            sumAngle -= LUT_CORDIC_ATAN[i];
        }

    }

    if (theta > (90 << CORDIC_MATH_FRACTION_BITS) &&
        theta < (270 << CORDIC_MATH_FRACTION_BITS)) {
        x = -x;    
        y = -y;
    }

    return (y << CORDIC_MATH_FRACTION_BITS) / x;
}

int32_t cordic_rectangular_polar(Coordinates *input) {
    int tempX, sumAngle = 0, x = input->x, y = input->y;
    if (x < 0 && y >= 0) {
        sumAngle = 90 * (1 << CORDIC_MATH_FRACTION_BITS);
        x = abs(x);
    } else if (x < 0 && y < 0) {
        sumAngle = 180 * (1 << CORDIC_MATH_FRACTION_BITS);
        x = abs(x);
        y = abs(y);
    }

    for (int i = 0; i < 15; i++) {
        tempX = x;
        if (y > 0) {
            /* Rotate clockwise */
            x += (y >> i);
            y -= (tempX >> i);
            sumAngle += LUT_CORDIC_ATAN[i];
        } else {
            /* Rotate counterclockwise */
            x -= (y >> i);
            y += (tempX >> i);
            sumAngle -= LUT_CORDIC_ATAN[i];
        }
    }
    input->theta = sumAngle;
    input->r = ((long)x * CORDIC_GAIN) >> CORDIC_MATH_FRACTION_BITS;
    return 0;
}

int32_t cordic_polar_rectangular(Coordinates *input) {
    int tempX, sumAngle = 0, x = CORDIC_GAIN, y = 0;

    input->theta %= (360 << CORDIC_MATH_FRACTION_BITS);

    if (input->theta > (90 * (1 << CORDIC_MATH_FRACTION_BITS)) &&
        input->theta < (270 * (1 << CORDIC_MATH_FRACTION_BITS))) {
        sumAngle = 180 * (1 << CORDIC_MATH_FRACTION_BITS);
        x = -x;
    }
    for (int i = 0; i < 15; i++) {
        tempX = x;
        if (input->theta > sumAngle) {
            /* Rotate counter clockwise */
            x -= (y >> i);
            y += (tempX >> i);
            sumAngle += LUT_CORDIC_ATAN[i];
        } else {
            /* Rotate clockwise */
            x += (y >> i);
            y -= (tempX >> i);
            sumAngle -= LUT_CORDIC_ATAN[i];
        }
    }
    input->x = ((long)x * input->r) >> CORDIC_MATH_FRACTION_BITS;
    input->y = ((long)y * input->r) >> CORDIC_MATH_FRACTION_BITS;
    return 0;
}

/*****************************Hyperbolic functions*****************************/

int32_t cordic_sqrt(int32_t x) {
    int poweroftwo;
    int y;

    if (x == 0) {
        return 0;
    }
    if (x == DECIMAL_TO_FP) {
        return DECIMAL_TO_FP;
    }
    poweroftwo = DECIMAL_TO_FP;

    if (x < DECIMAL_TO_FP) {
        while (x <= (((long)poweroftwo * poweroftwo) >> CORDIC_MATH_FRACTION_BITS)) {
            poweroftwo >>= 1;
        }
        y = poweroftwo;
    } else if (x > DECIMAL_TO_FP) {
        while ((((long)poweroftwo * poweroftwo) >> CORDIC_MATH_FRACTION_BITS) <= x) {
            poweroftwo <<= 1;
        }
        y = poweroftwo >> 1;
    }
    for (int i = 1; i <= 15; i++) {
        poweroftwo >>= 1;
        if (((long)(y + poweroftwo) * (y + poweroftwo) >> CORDIC_MATH_FRACTION_BITS) <= x) {
            y = y + poweroftwo;
        }
    }
    return y;
}

int32_t cordic_arctanh(int32_t y, int32_t x) {
    int tempX, k = 4, sumAngle = 0;

    for (int i = 1; i < 15; i++) {
        tempX = x;
        if (y < 0) {
            /* Rotate clockwise */
            x += (y >> i);
            y += (tempX >> i);
            sumAngle -= LUT_CORDIC_ATANH[i - 1];
        } else {
            /* Rotate counterclockwise */
            x -= (y >> i);
            y -= (tempX >> i);
            sumAngle += LUT_CORDIC_ATANH[i - 1];
        }
        if (i == k) {
            k = (3 * k) + 1;
            tempX = x;
            if (y < 0) {
                /* Rotate clockwise */
                x += (y >> i);
                y += (tempX >> i);
                sumAngle -= LUT_CORDIC_ATANH[i - 1];
            } else {
                /* Rotate counterclockwise */
                x -= (y >> i);
                y -= (tempX >> i);
                sumAngle += LUT_CORDIC_ATANH[i - 1];
            }
        }
    }
    return sumAngle;
}

int32_t cordic_ln(int32_t input) {
    int k = 0;
    long calculate = input;

    while (calculate > EULER) {
        calculate <<= CORDIC_MATH_FRACTION_BITS;
        //printf("calculate : %ld\n", calculate);
        calculate /= EULER;
        k += DECIMAL_TO_FP;
    }

    int y = calculate - DECIMAL_TO_FP;
    int x = calculate + DECIMAL_TO_FP;

    return (to_radians(cordic_arctanh(y,x) << 1) + k);
}

int32_t cordic_arccosh(int32_t x) {
    int tempX, k = 4, sumAngle = 0, y = DECIMAL_TO_FP, xt = x;

    for (int i = 1; i < 15; i++) {
        tempX = x;
        if (y < 0) {
            /* Rotate clockwise */
            x += (y >> i);
            y += (tempX >> i);
            sumAngle -= LUT_CORDIC_ATANH[i - 1];
        } else {
            /* Rotate counterclockwise */
            x -= (y >> i);
            y -= (tempX >> i);
            sumAngle += LUT_CORDIC_ATANH[i - 1];
        }
        if (i == k) {
            k = (3 * k) + 1;
            tempX = x;
            if (y < 0) {
                /* Rotate clockwise */
                x += (y >> i);
                y += (tempX >> i);
                sumAngle -= LUT_CORDIC_ATANH[i - 1];
            } else {
                /* Rotate counterclockwise */
                x -= (y >> i);
                y -= (tempX >> i);
                sumAngle += LUT_CORDIC_ATANH[i - 1];
            }
        }
    }

    return to_degree(cordic_ln((( (long)x << CORDIC_MATH_FRACTION_BITS) / CORDIC_GAIN_HYPERBOLIC_VECTOR) + xt));
}

int32_t cordic_arcsinh(int32_t y) {
    int tempX, k = 4, sumAngle = 0, x = DECIMAL_TO_FP, yt = y;

    for (int i = 0; i < 15; i++) {
        tempX = x;
        if (y < 0) {
            /* Rotate clockwise */
            x -= (y >> i);
            y += (tempX >> i);
            sumAngle -= LUT_CORDIC_ATAN[i];
        } else {
            /* Rotate counterclockwise */
            x += (y >> i);
            y -= (tempX >> i);
            sumAngle += LUT_CORDIC_ATAN[i];
        }
    }

    return to_degree(cordic_ln((((long)x << CORDIC_MATH_FRACTION_BITS) / CORDIC_GAIN_HYPERBOLIC_CIRCULAR) + yt));
}

int32_t cordic_sinh(int32_t theta) {
    int tempX, k = 4, sumAngle = theta, y = 0,
               x = ONE_DIV_CORDIC_GAIN_HYPERBOLIC;
    for (int i = 1; i < 15; i++) {
        tempX = x;
        if (sumAngle > 0) {
            /* Rotate clockwise */
            x += (y >> i);
            y += (tempX >> i);
            sumAngle -= LUT_CORDIC_ATANH[i - 1];
        } else {
            /* Rotate counterclockwise */
            x -= (y >> i);
            y -= (tempX >> i);
            sumAngle += LUT_CORDIC_ATANH[i - 1];
        }
        if (i == k) {
            k = (3 * k) + 1;
            tempX = x;
            if (sumAngle > 0) {
                /* Rotate clockwise */
                x += (y >> i);
                y += (tempX >> i);
                sumAngle -= LUT_CORDIC_ATANH[i - 1];
            } else {
                /* Rotate counterclockwise */
                x -= (y >> i);
                y -= (tempX >> i);
                sumAngle += LUT_CORDIC_ATANH[i - 1];
            }
        }
    }
    return y;
}

int32_t cordic_cosh(int32_t theta) {
    int tempX, k = 4, sumAngle = theta, y = 0,
               x = ONE_DIV_CORDIC_GAIN_HYPERBOLIC;
    for (int i = 1; i < 15; i++) {
        tempX = x;
        if (sumAngle > 0) {
            /* Rotate clockwise */
            x += (y >> i);
            y += (tempX >> i);
            sumAngle -= LUT_CORDIC_ATANH[i - 1];
        } else {
            /* Rotate counterclockwise */
            x -= (y >> i);
            y -= (tempX >> i);
            sumAngle += LUT_CORDIC_ATANH[i - 1];
        }
        if (i == k) {
            k = (3 * k) + 1;
            tempX = x;
            if (sumAngle > 0) {
                /* Rotate clockwise */
                x += (y >> i);
                y += (tempX >> i);
                sumAngle -= LUT_CORDIC_ATANH[i - 1];
            } else {
                /* Rotate counterclockwise */
                x -= (y >> i);
                y -= (tempX >> i);
                sumAngle += LUT_CORDIC_ATANH[i - 1];
            }
        }
    }
    return x;
}

int32_t cordic_tanh(int32_t theta) {
    return (((long)cordic_sinh(theta) << CORDIC_MATH_FRACTION_BITS) / cordic_cosh(theta));
}

int32_t cordic_exp(int32_t exponent) {
    int tempX, k = 4, sumAngle = to_degree(exponent),
               y = ONE_DIV_CORDIC_GAIN_HYPERBOLIC,
               x = ONE_DIV_CORDIC_GAIN_HYPERBOLIC, n = 0;

    while (sumAngle > ONE_EIGHTY_DIV_PI) {
        sumAngle -= ONE_EIGHTY_DIV_PI;
        n++;
    }

    for (int i = 1; i < 15; i++) {
        tempX = x;
        if (sumAngle > 0) {
            /* Rotate clockwise */
            x += (y >> i);
            y += (tempX >> i);
            sumAngle -= LUT_CORDIC_ATANH[i - 1];
        } else {
            /* Rotate counterclockwise */
            x -= (y >> i);
            y -= (tempX >> i);
            sumAngle += LUT_CORDIC_ATANH[i - 1];
        }
        if (i == k) {
            k = (3 * k) + 1;
            tempX = x;
            if (sumAngle > 0) {
                /* Rotate clockwise */
                x += (y >> i);
                y += (tempX >> i);
                sumAngle -= LUT_CORDIC_ATANH[i - 1];
            } else {
                /* Rotate counterclockwise */
                x -= (y >> i);
                y -= (tempX >> i);
                sumAngle += LUT_CORDIC_ATANH[i - 1];
            }
        }
    }

    y = DECIMAL_TO_FP;
    for (int i = 0; i < n; i++) {
        y = ((long)y * EULER) >> CORDIC_MATH_FRACTION_BITS;
    }
    return ((long)x * y) >> CORDIC_MATH_FRACTION_BITS;
}

int32_t cordic_pow(int32_t base, int32_t exponent) {
    return cordic_exp(exponent * cordic_ln(base) >> CORDIC_MATH_FRACTION_BITS);
}

/*****************************None Cordic Math*****************************/

int32_t abs(int32_t input) {
    if (input >= 0) {
        return input;
    } else {
        return -input;
    }
}

int32_t isEven(int32_t input) {
    if (input % 2) return 0;
    return 1;
}

int32_t isOdd(int32_t input) {
    if (input % 2) return 1;
    return 0;
}

int32_t to_degree(int32_t input) {
    return ((long)input * ONE_EIGHTY_DIV_PI >> CORDIC_MATH_FRACTION_BITS);
}

int32_t to_radians(int32_t input) {
    return (((long)input << CORDIC_MATH_FRACTION_BITS) / ONE_EIGHTY_DIV_PI);
}
