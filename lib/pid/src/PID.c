#include "PID.h"

void PIDController_Init(PIDController *pid) {

	pid->Kp = 3.50f; //3.5
    pid->Kd = 0.1f; //40.0
    pid->Ki = 40.0f; //0.10
	/* Clear controller variables */
	pid->integrator = 0.0f;
	pid->prevError  = 0.0f;

	pid->limMax = 230400;
	pid->limMin = -230400;

	pid->limMaxInt = 100;
	pid->limMinInt = -100;
	

	pid->differentiator  = 0.0f;
	pid->prevMeasurement = 0.0f;

	pid->out = 0.0f;

}

int millis(void){
	uint64_t mtime = get_timer_value();
	return ((mtime*4000.0)/SystemCoreClock);
}

float PIDController_Update(PIDController *pid, float setpoint, float measurement) { 
	/*
	* Error signal
	*/
    float error = setpoint - measurement; //setpoint = börvärde, measurement = ärvärde

	/*
	* Proportional
	*/
    float proportional = pid->Kp * error;

	/*
	* Integral
	*/
    pid->integrator += (pid->Ki * error);
	
	/* Anti-wind-up via integrator clamping */
    if (pid->integrator > pid->limMaxInt) {

        pid->integrator = pid->limMaxInt;

    } else if (pid->integrator < pid->limMinInt) {

        pid->integrator = pid->limMinInt;

    }

	/*
	* Derivative (band-limited differentiator)
	*/
		
    pid->differentiator = -(pid->Kd * (measurement - pid->prevMeasurement));
	/*
	* Compute output and apply limits
	*/
    pid->out = proportional + pid->integrator + pid->differentiator;

    if (pid->out > pid->limMax) {

        pid->out = pid->limMax;

    } else if (pid->out < pid->limMin) {

        pid->out = pid->limMin;

    }

	/* Store error and measurement for later use */
    pid->prevError       = error;
    pid->prevMeasurement = measurement;

	/* Return controller output */
    return pid->out;

}