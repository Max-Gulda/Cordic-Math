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
#include "cordic-math.h"

/*****************************************VECTORING MODE***********************************************/
int32_t cordic_atan(int32_t y, int32_t x){
    int sumAngle = 0, tempX;
    if(x<0){
        x = -x;
        y = -y;
    }
    for (int i = 0; i < 15; i++){
        tempX = x;
        if(y>0){
            /* Rotate clockwise */
            x += (y>>i);
            y -= (tempX>>i);
            sumAngle += LUT_CORDIC_ATAN[i];
        }else{
            /* Rotate counterclockwise */
            x -= (y>>i);
            y += (tempX>>i);
            sumAngle -= LUT_CORDIC_ATAN[i];
        }
    }
    return sumAngle;
}


int32_t cordic_hypotenuse(int32_t y, int32_t x){
    int tempX;
    x = abs(x);
    y = abs(y);

    for (int i = 0; i < 15; i++){
        tempX = x;
        if(y>0){
            /* Rotate clockwise */
            x += (y>>i);
            y -= (tempX>>i);
        }else{
            /* Rotate counterclockwise */
            x -= (y>>i);
            y += (tempX>>i);
        }
    }
    return (x*CORDIC_GAIN)>>8;
}


int32_t cordic_cos(int32_t theta){
    int x = CORDIC_GAIN, y = 0, sumAngle = 0, tempX;

    theta %= (360*DECIMAL_TO_FP);
    if(theta>(90*DECIMAL_TO_FP)){
        sumAngle = 180*DECIMAL_TO_FP;
    }
    if(theta>(270*DECIMAL_TO_FP)){
        sumAngle = 360*DECIMAL_TO_FP;
    }

    for(int i = 0; i < 15; i++){
        tempX = x;
        if(theta > sumAngle){
            /* Rotate counter clockwise */
            x -= (y>>i);
            y += (tempX>>i);
            sumAngle += LUT_CORDIC_ATAN[i];
        }else{
            /* Rotate clockwise */
            x += (y>>i);
            y -= (tempX>>i);
            sumAngle -= LUT_CORDIC_ATAN[i];
        }
    }
    if(theta > (90*DECIMAL_TO_FP) && theta < (270*DECIMAL_TO_FP)){
        x = -x;
    }
    return x;
}


int32_t cordic_sin(int32_t theta){
    int x = CORDIC_GAIN, y = 0, sumAngle = 0, tempX;

    theta %= (360*DECIMAL_TO_FP);
    
    if(theta>(90*DECIMAL_TO_FP)){
        sumAngle = 180*DECIMAL_TO_FP;
    }
    if(theta>(270*DECIMAL_TO_FP)){
        sumAngle = 360*DECIMAL_TO_FP;
    }

    for(int i = 0; i < 15; i++){
        tempX = x;
        if(theta > sumAngle){
            /* Rotate counter clockwise */
            x -= (y>>i);
            y += (tempX>>i);
            sumAngle += LUT_CORDIC_ATAN[i];
        }else{
            /* Rotate clockwise */
            x += (y>>i);
            y -= (tempX>>i);
            sumAngle -= LUT_CORDIC_ATAN[i];
        }
    }
    if(theta < (180*DECIMAL_TO_FP) || theta > (360*DECIMAL_TO_FP)){
        y = abs(y);
    }
    return y;
}


int32_t cordic_asin(int32_t yInput){
    int x = CORDIC_GAIN, y = 0, sumAngle = 0, tempX;

    for(int i = 0; i < 15; i++){
        tempX = x;
        if(y < yInput){
            /* Rotate counter clockwise */
            x -= (y>>i);
            y += (tempX>>i);
            sumAngle += LUT_CORDIC_ATAN[i];
        }else{
            /* Rotate clockwise */
            x += (y>>i);
            y -= (tempX>>i);
            sumAngle -= LUT_CORDIC_ATAN[i];
        }
    }
    if(sumAngle<-90*DECIMAL_TO_FP){
        sumAngle=-90*DECIMAL_TO_FP;
    }else if(sumAngle>90*DECIMAL_TO_FP){
        sumAngle=90*DECIMAL_TO_FP;
    }
    return sumAngle;
}


int32_t cordic_acos(int32_t xInput){
    int x = 0, y = CORDIC_GAIN, sumAngle = 90*DECIMAL_TO_FP, tempX;

    for(int i = 0; i < 15; i++){
        tempX = x;
        if(x > xInput){
            /* Rotate counter clockwise */
            x -= (y>>i);
            y += (tempX>>i);
            sumAngle += LUT_CORDIC_ATAN[i];
        }else{
            /* Rotate clockwise */
            x += (y>>i);
            y -= (tempX>>i);
            sumAngle -= LUT_CORDIC_ATAN[i];
        }
    }
    if(sumAngle>180*DECIMAL_TO_FP){
        sumAngle = 180*DECIMAL_TO_FP;
    }else if(sumAngle<0){
        sumAngle=0;
    }
    return sumAngle;
}


int32_t cordic_tan(int32_t degree){
    return (cordic_sin(degree)*DECIMAL_TO_FP/cordic_cos(degree));
}

/*****************************************HYPERBOLIC MODE***********************************************/

/**
 * @brief Fast fixedpoint [24|8] squareroot using the cordic algorithm
 * 
 * @param x, 
 * @return int32_t 
 */

int32_t cordic_sqrt(int32_t x){
    int poweroftwo;
    int y;

    if(x==0){
        return 0;
    }
    if(x==DECIMAL_TO_FP){
        return DECIMAL_TO_FP;
    }
    poweroftwo = DECIMAL_TO_FP;

    if(x < DECIMAL_TO_FP){
        while(x <= ((poweroftwo * poweroftwo) >> 8)){
            poweroftwo >>= 1;
        }
        y = poweroftwo;
    }else if(x > DECIMAL_TO_FP){
        while(((poweroftwo * poweroftwo) >> 8) <= x){
            poweroftwo <<= 1;
        }
        y = poweroftwo >> 1;
    }
    for(int i = 1; i <= 15; i++){
        poweroftwo >>= 1;
        if(((y + poweroftwo ) * ( y + poweroftwo ) >> 8) <= x){
            y = y + poweroftwo;
        }
    }
    return y;
}


int32_t cordic_arctanh(int32_t y, int32_t x){
    int tempX, k = 4, sumAngle=0;
    
    for (int i = 1; i < 15; i++){
        tempX = x;
        if(y<0){
            /* Rotate clockwise */
            x += (y>>i);
            y += (tempX>>i);
            sumAngle -= LUT_CORDIC_ATANH[i-1];
        }else{
            /* Rotate counterclockwise */
            x -= (y>>i);
            y -= (tempX>>i);
            sumAngle += LUT_CORDIC_ATANH[i-1];
        }
        if(i==k){
            k = (3*k) + 1;
            tempX = x;
            if(y<0){
                /* Rotate clockwise */
                x += (y>>i);
                y += (tempX>>i);
                sumAngle -= LUT_CORDIC_ATANH[i-1];
            }else{
                /* Rotate counterclockwise */
                x -= (y>>i);
                y -= (tempX>>i);
                sumAngle += LUT_CORDIC_ATANH[i-1];
            }
        }
    }
    return sumAngle;
}


int32_t cordic_ln(int32_t input){
    int k = 0;
    while(input>EULER){
        input<<=8;
        input/=EULER;
        k+=DECIMAL_TO_FP;
    }
    return (to_radians(cordic_arctanh(input-DECIMAL_TO_FP,input+DECIMAL_TO_FP))<<1)+k;
}


int32_t cordic_arccosh(int32_t x){
    int tempX, k = 4, sumAngle=0, y = DECIMAL_TO_FP, xt = x;
    
    for (int i = 1; i < 15; i++){
        tempX = x;
        if(y<0){
            /* Rotate clockwise */
            x += (y>>i);
            y += (tempX>>i);
            sumAngle -= LUT_CORDIC_ATANH[i-1];
        }else{
            /* Rotate counterclockwise */
            x -= (y>>i);
            y -= (tempX>>i);
            sumAngle += LUT_CORDIC_ATANH[i-1];
        }
        if(i==k){
            k = (3*k) + 1;
            tempX = x;
            if(y<0){
                /* Rotate clockwise */
                x += (y>>i);
                y += (tempX>>i);
                sumAngle -= LUT_CORDIC_ATANH[i-1];
            }else{
                /* Rotate counterclockwise */
                x -= (y>>i);
                y -= (tempX>>i);
                sumAngle += LUT_CORDIC_ATANH[i-1];
            }
        }
    }
    
    return to_degree(cordic_ln(((x<<8)/212) + xt));
}


int32_t cordic_arcsinh(int32_t y){
    int tempX, k = 4, sumAngle=0, x = DECIMAL_TO_FP, yt = y;
    
    for (int i = 0; i < 15; i++){
        tempX = x;
        if(y<0){
            /* Rotate clockwise */
            x -= (y>>i);
            y += (tempX>>i);
            sumAngle -= LUT_CORDIC_ATAN[i];
        }else{
            /* Rotate counterclockwise */
            x += (y>>i);
            y -= (tempX>>i);
            sumAngle += LUT_CORDIC_ATAN[i];
        }
    }
    
    return to_degree(cordic_ln(((x << 8)/422) + yt));
}


int32_t cordic_sinh(int32_t theta){
    int tempX, k = 4, sumAngle = theta, y = 0, x = 309;
    for (int i = 1; i < 15; i++){
        tempX = x;
        if(sumAngle > 0){
            /* Rotate clockwise */
            x += (y>>i);
            y += (tempX>>i);
            sumAngle -= LUT_CORDIC_ATANH[i-1];
        }else{
            /* Rotate counterclockwise */
            x -= (y>>i);
            y -= (tempX>>i);
            sumAngle += LUT_CORDIC_ATANH[i-1];
        }
        if(i==k){
            k = (3*k) + 1;
            tempX = x;
            if(sumAngle > 0){
                /* Rotate clockwise */
                x += (y>>i);
                y += (tempX>>i);
                sumAngle -= LUT_CORDIC_ATANH[i-1];
            }else{
                /* Rotate counterclockwise */
                x -= (y>>i);
                y -= (tempX>>i);
                sumAngle += LUT_CORDIC_ATANH[i-1];
            }
        }
    }
    return y;
}


int32_t cordic_cosh(int32_t theta){
    int tempX, k = 4, sumAngle = theta, y = 0, x = 309;
    for (int i = 1; i < 15; i++){
        tempX = x;
        if(sumAngle > 0){
            /* Rotate clockwise */
            x += (y>>i);
            y += (tempX>>i);
            sumAngle -= LUT_CORDIC_ATANH[i-1];
        }else{
            /* Rotate counterclockwise */
            x -= (y>>i);
            y -= (tempX>>i);
            sumAngle += LUT_CORDIC_ATANH[i-1];
        }
        if(i==k){
            k = (3*k) + 1;
            tempX = x;
            if(sumAngle > 0){
                /* Rotate clockwise */
                x += (y>>i);
                y += (tempX>>i);
                sumAngle -= LUT_CORDIC_ATANH[i-1];
            }else{
                /* Rotate counterclockwise */
                x -= (y>>i);
                y -= (tempX>>i);
                sumAngle += LUT_CORDIC_ATANH[i-1];
            }
        }
    }
    return x;
}


int32_t cordic_tanh(int32_t theta){
    return ((cordic_sinh(theta) << 8)/ cordic_cosh(theta));
}


int32_t cordic_exp(int32_t exponent){
    int tempX, k = 4, sumAngle = to_degree(exponent), y = 309, x = 309, n = 0;

    while(sumAngle>14667){
        sumAngle-=14667;
        n++;
    }

    for (int i = 1; i < 15; i++){
        tempX = x;
        if(sumAngle > 0){
            /* Rotate clockwise */
            x += (y>>i);
            y += (tempX>>i);
            sumAngle -= LUT_CORDIC_ATANH[i-1];
        }else{
            /* Rotate counterclockwise */
            x -= (y>>i);
            y -= (tempX>>i);
            sumAngle += LUT_CORDIC_ATANH[i-1];
        }
        if(i==k){
            k = (3*k) + 1;
            tempX = x;
            if(sumAngle > 0){
                /* Rotate clockwise */
                x += (y>>i);
                y += (tempX>>i);
                sumAngle -= LUT_CORDIC_ATANH[i-1];
            }else{
                /* Rotate counterclockwise */
                x -= (y>>i);
                y -= (tempX>>i);
                sumAngle += LUT_CORDIC_ATANH[i-1];
            }
        }
    }
    y = DECIMAL_TO_FP;
    for(int i = 0; i < n; i++){
        y = (y * EULER) >> 8;
    }
    return (x*y)>>8;
}


int32_t cordic_pow(int32_t base, int32_t exponent){
    return cordic_exp(exponent*cordic_ln(base) >> 8);
}


int32_t abs(int32_t input){
    if(input>0){
        return input;
    }else{
        return -input;
    }
}


int32_t isEven(int32_t input){
    if(input%2) return 0;
    return 1;
}


int32_t isOdd(int32_t input){
    if(input%2) return 1;
    return 0;
}


int32_t to_degree(int32_t input){
    return (input*14667>>8);
}


int32_t to_radians(int32_t input){
    return ((input<<8)/14667);
}