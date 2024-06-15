#include <stdio.h>
#include <stdint.h>

#define SCALE (1 << 8)
#define CORDIC_MATH_FRACTION_BITS 8
#define CORDIC_SPEED_FACTOR 15

// Fixed-point constants
#define FLOAT_TO_INT(x) ((x) >= 0 ? (int)((x) + 0.5) : (int)((x)-0.5))
#define CORDIC_GAIN FLOAT_TO_INT(0.607252935 * SCALE)

// Look-Up Table for arctan values
static const int32_t LUT_CORDIC_ATAN[CORDIC_SPEED_FACTOR] = {
    FLOAT_TO_INT(45.0 * (1 << CORDIC_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(26.5651 * (1 << CORDIC_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(14.0362 * (1 << CORDIC_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(7.1250 * (1 << CORDIC_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(3.5763 * (1 << CORDIC_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(1.7899 * (1 << CORDIC_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(0.8952 * (1 << CORDIC_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(0.4476 * (1 << CORDIC_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(0.2238 * (1 << CORDIC_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(0.1119 * (1 << CORDIC_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(0.0560 * (1 << CORDIC_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(0.0280 * (1 << CORDIC_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(0.0140 * (1 << CORDIC_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(0.0070 * (1 << CORDIC_MATH_FRACTION_BITS)),
    FLOAT_TO_INT(0.0035 * (1 << CORDIC_MATH_FRACTION_BITS))
};

int32_t cordic_abs(int32_t input) {
    return (input >= 0) ? input : -input;
}

int32_t cordic_atan(int32_t y, int32_t x) {
    int sumAngle = 0, tempX;
    if (x < 0) {
        x = -x;
        y = -y;
    }
    for (int i = 0; i < CORDIC_SPEED_FACTOR; i++) {
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

/**
 * @brief Fast fixed-point calculation of hypotenuse using the CORDIC algorithm
 *
 * @param y fixed-point according to CORDIC_MATH_FRACTION_BITS
 * @param x fixed-point according to CORDIC_MATH_FRACTION_BITS
 *
 * @return 32-bit int fixed-point according to CORDIC_MATH_FRACTION_BITS, sqrt(x*x + y*y)
 */
int32_t cordic_hypotenuse(int32_t y, int32_t x) {
    int tempX;
    x = cordic_abs(x);
    y = cordic_abs(y);

    for (int i = 0; i < CORDIC_SPEED_FACTOR; i++) {
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

    return ((int64_t)x * CORDIC_GAIN) >> CORDIC_MATH_FRACTION_BITS;
}

// Fixed-point multiplication
int32_t fixed_mul(int32_t a, int32_t b) {
    return (a * b) / SCALE;
}

// Fixed-point division
int32_t fixed_div(int32_t a, int32_t b) {
    return (a * SCALE) / b;
}

// Structure to hold Kalman filter state
typedef struct {
    int32_t q; // Process noise covariance
    int32_t r; // Measurement noise covariance
    int32_t x; // Estimated value
    int32_t p; // Estimation error covariance
    int32_t k; // Kalman gain
} KalmanFilter;

// Initialize the Kalman filter
void kalman_init(KalmanFilter *kf, int32_t q, int32_t r, int32_t initial_value) {
    kf->q = q;
    kf->r = r;
    kf->x = initial_value;
    kf->p = SCALE; // Initialize p with 1 in fixed-point
    kf->k = 0;
}

// Update the Kalman filter with a new measurement
void kalman_update(KalmanFilter *kf, int32_t measurement) {
    // Prediction update
    kf->p = kf->p + kf->q;

    // Measurement update
    kf->k = fixed_div(kf->p, kf->p + kf->r);
    kf->x = kf->x + fixed_mul(kf->k, measurement - kf->x);
    kf->p = fixed_mul((SCALE - kf->k), kf->p);
}

// Function to calculate pitch and roll from accelerometer data using CORDIC
void calculate_angles(int32_t ax, int32_t ay, int32_t az, int32_t *pitch, int32_t *roll) {
    int32_t hyp = cordic_hypotenuse(ay, az);
    *pitch = cordic_atan(-ax, hyp);
    *roll = cordic_atan(ay, az);
}

// Example usage
int main() {
    // Example accelerometer measurements for multiple runs (scaled)
    int32_t measurements[][3] = {
        {0.0 * SCALE, 0.0 * SCALE, 1.0 * SCALE},
        {0.1 * SCALE, 0.1 * SCALE, 0.98 * SCALE},
        {0.2 * SCALE, 0.1 * SCALE, 0.97 * SCALE},
        {-0.1 * SCALE, 0.1 * SCALE, 0.99 * SCALE},
        {-0.2 * SCALE, -0.1 * SCALE, 0.98 * SCALE},
        {0.0 * SCALE, -0.2 * SCALE, 0.96 * SCALE}
    };
    int num_measurements = sizeof(measurements) / sizeof(measurements[0]);

    // Initialize the Kalman filters for pitch and roll
    KalmanFilter kf_pitch;
    KalmanFilter kf_roll;
    kalman_init(&kf_pitch, 0.1 * SCALE, 0.1 * SCALE, 0.0 * SCALE);
    kalman_init(&kf_roll, 0.1 * SCALE, 0.1 * SCALE, 0.0 * SCALE);

    // Process each set of measurements
    for (int i = 0; i < num_measurements; i++) {
        int32_t ax = measurements[i][0];
        int32_t ay = measurements[i][1];
        int32_t az = measurements[i][2];

        // Calculate pitch and roll from accelerometer
        int32_t pitch, roll;
        calculate_angles(ax, ay, az, &pitch, &roll);

        // Update Kalman filters with the calculated pitch and roll
        kalman_update(&kf_pitch, pitch);
        kalman_update(&kf_roll, roll);

        // Print results
        printf("Run %d\n", i + 1);
        printf("Measurement: ax=%f, ay=%f, az=%f\n", (double)ax / SCALE, (double)ay / SCALE, (double)az / SCALE);
        printf("Pitch: %f, Estimated Pitch: %f\n", (double)pitch / SCALE, (double)kf_pitch.x / SCALE);
        printf("Roll: %f, Estimated Roll: %f\n", (double)roll / SCALE, (double)kf_roll.x / SCALE);
        printf("\n");
    }

    return 0;
}
