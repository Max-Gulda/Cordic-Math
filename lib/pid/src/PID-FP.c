#include "PID-FP.h"

void PIDControllerFP_Init(PIDController *pid) {

	/*
	* Kp, Ki and Kd is Kx * 256
	*/

	pid->Kp = 896;
    pid->Ki = 10240;
	pid->Kd = 25;
	/* Clear controller variables */
	pid->proportional = 0;
	
	pid->integrator = 0;
	pid->prevError  = 0;

	pid->limMax = 230400; //900 << 8
	pid->limMin = -230400;

	pid->limMaxInt = 2949120; //45 << 16
	pid->limMinInt = -2949120;

	pid->differentiator  = 0;
	pid->prevMeasurement = 0;

	pid->out = 0;

}

int millis(void){
	uint64_t mtime = get_timer_value();
	return ((mtime*4000.0)/SystemCoreClock);
}
/**
 * @brief PID-controller
 * 
 * @param pid PID struct
 * @param setpoint Setpoint fixedpoint << 8
 * @param measurement measurement fixedpoint << 8
 * @param deltaTime time since last loop in mS
 * @return float 
 */
int PIDControllerFP_Update(PIDController *pid, int setpoint, int measurement, int deltaTime) { 
	/*
	* Error signal
	*/
    int error = setpoint - measurement; //setpoint = börvärde, measurement = ärvärde

	/*
	* Proportional
	*/
    pid->proportional = pid->Kp * error;
	/*
	* Integral
	*/
    pid->integrator += (pid->Ki * error*deltaTime/1000);
	
	/* Anti-wind-up via integrator clamping */
    if (pid->integrator > pid->limMaxInt) {

        pid->integrator = pid->limMaxInt;

    } else if (pid->integrator < pid->limMinInt) {

        pid->integrator = pid->limMinInt;

    }

	/*
	* Derivative (band-limited differentiator)
	*/
		
    pid->differentiator = -(pid->Kd * ((measurement - pid->prevMeasurement)*1000/deltaTime));

	/*
	* Compute output and apply limits
	*/
    pid->out = pid->proportional + pid->integrator + pid->differentiator;

    if (pid->out > pid->limMax) {

        pid->out = pid->limMax;

    } else if (pid->out < pid->limMin) {

        pid->out = pid->limMin;

    }
	pid->out = pid->out >> 8; /* 16/16 fixedpoint converted to 24/8 */

	/* Store error and measurement for later use */
    pid->prevError       = error;
    pid->prevMeasurement = measurement;

	/* Return controller output */
    return pid->out;

}
