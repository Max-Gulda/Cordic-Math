<h1 align="center">Cordic math</h3>

  <p align="center">
    Project Decription
    <br />
    <a href="https://github.com/Flaxyson/Cordic-Math"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    ·
    <a href="https://github.com/Flaxyson/Cordic-Math/issues">Report Bug</a>
    ·
    <a href="https://github.com/Flaxyson/Cordic-Math/issues">Request Feature</a>
  </p>
</div>

<br />

## About The CORDIC Algorithm

The cordic algorithm is a simple and efficient algorithm to calculate trigonometric functions, hyperbolic functions, square roots, multiplications, divisions, exponentials, logarithms and more. The cordic algorithm dates back to 1956 but still work great for microcontrollers and such today.

Cordic uses simple simple bit-shift operations for several computing tasks. In the library I included other mathimatical functions aswell that might be useful for future users.



<p align="right">(<a href="#top">back to top</a>)</p>


## Modifying The Code

The code is built around the defined variable CORDIC_MATH_FRACTION_BITS in cordic-math.h file. As default this is set to 16, which means the calculations is made in 16 bit fixed point arithmetic. This variable can easily be changed by future users to change the fixedpoint to suite their needs.

![my image](pictures/Screenshoot-1.png) 


## Accuracy
The following image shows the precision of the Cordic algorithm. The left column is math.h and the right column is the Cordic algorithm. The input was randomly picked and the result was the following:

![my image](pictures/Screenshoot-2.png)



## Functions

- [x] Arctan
- [x] Arcsin
- [x] Arccos
- [x] Tan
- [x] Cos
- [x] Sin
- [x] Squareroot
- [x] Calculation of Hypotenuse
- [x] Arctan Hyperbolic
- [x] Arcsin Hyperbolic
- [x] Arccos Hyperbolic
- [x] Tan Hyperbolic
- [x] Cos Hyperbolic
- [x] Sin Hyperbolic
- [x] e to the Power
- [x] x to the Power
- [x] Natural logratihm
- [x] Conversion From Radians to Degrees
- [x] Conversion From Degrees to Radians
- [x] Absolute
- [x] Is Odd
- [x] Is Even
- [ ] Rectangular to Polar Conversion
- [ ] Polar to Rectangular Conversion


<p align="right">(<a href="#top">back to top</a>)</p>

